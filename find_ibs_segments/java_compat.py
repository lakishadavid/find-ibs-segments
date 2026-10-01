# find-ibs-segments.  This module comes from merge-ibs-segments, a modified
# Python translation of merge-ibd-segments from Refined-IBD version 17Jan20.102.
#
# Copyright (C) 2014-2016 Brian L. Browning
# Copyright (C) 2026 LaKisha T. David
#
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
# more details.
#
# You should have received a copy of the GNU General Public License along
# with this program.  If not, see <https://www.gnu.org/licenses/>.
#
# Modified by LaKisha T. David, September 29-30, 2026: translated from Java to
# Python from src/blbutil/StringUtil.java and src/blbutil/Utilities.java in
# refined-ibd.17Jan20.102.zip (SHA-256
# 0b0abf48528ec53d3fec7f21a9742c6e5960857c5a225b0e49346179bc45e6ea). The Java
# SE library behaviour in this file is not translated from any source code; see
# the module documentation.  Later changes are recorded in the Git history.
#
# Copied into find-ibs-segments on October 1, 2026, with only the first lines
# of this notice and this paragraph changed, from
# merge-ibs-segments (https://github.com/lakishadavid/merge-ibs-segments),
# commit 93402c9.  Later changes are recorded in the Git history.
#
# Brian L. Browning's original notice, reproduced unchanged from those files:
#
# Copyright (C) 2014-2016 Brian L. Browning
#
# This file is part of Refined-IBD
#
# Refined-IBD is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Beagle is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""Java behaviour that merge-ibd-segments depends on.

Two kinds of code are in this module.

- Translated from Brian Browning's Java: the field splitting of
  blbutil.StringUtil and the ways blbutil.Utilities ends the program.
- Java SE library behaviour that the jar shows: string trimming, character
  classes, Integer.parseInt, Double.parseDouble, Float.parseFloat,
  Arrays.binarySearch, Float.compare and DecimalFormat("0.##") /
  ("0.###").  These functions are written from the Java SE API
  documentation (Javadoc) and from the jar's observed output, recorded by
  the exact-conversion test in tests/parity (Java 25.0.4.1).  They are not
  translated from OpenJDK source code, whose licence (GPL version 2 only,
  with the Classpath Exception) does not allow its combination with this
  GPL version 3 program.

Strings are handled as Java handles them: as sequences of UTF-16 code units.
``to_java_chars`` splits every character outside the Basic Multilingual Plane
into its two surrogate code units, so that lengths, indices and per-character
tests agree with Java; ``java_encode`` joins them again for output.
"""

from __future__ import annotations

import math
import re
import struct
import unicodedata
from decimal import ROUND_HALF_EVEN, Decimal, localcontext
from fractions import Fraction
from typing import List, Optional, Tuple

INT_MIN = -(1 << 31)
INT_MAX = (1 << 31) - 1
DOUBLE_MAX = 1.7976931348623157e308


# ---------------------------------------------------------------------------
# How the Java program ends (blbutil.Utilities and uncaught exceptions)
# ---------------------------------------------------------------------------

class JavaException(Exception):
    """A Java exception.

    When one reaches ``main`` uncaught, the Java virtual machine writes
    ``Exception in thread "main" `` followed by the exception's class name,
    ``: `` and its message, then a stack trace, to standard error and exits
    with status 1.
    """

    def __init__(self, java_class: str, message: Optional[str] = None):
        super().__init__(java_class if message is None
                         else java_class + ": " + message)
        self.java_class = java_class
        self.message = message

    def to_string(self) -> str:
        """The class name, and ``: `` and the message if there is one."""
        if self.message is None:
            return self.java_class
        return self.java_class + ": " + self.message


class JavaTermination(Exception):
    """A deliberate ``System.exit``: text for standard output and standard
    error, and the exit status."""

    def __init__(self, stdout: str = "", stderr: str = "", status: int = 1):
        super().__init__(stderr or stdout)
        self.stdout = stdout
        self.stderr = stderr
        self.status = status


def illegal_argument(message: Optional[str]) -> JavaException:
    return JavaException("java.lang.IllegalArgumentException", message)


def number_format_for_input(s: str) -> JavaException:
    """The NumberFormatException the jar reports for a string that is not a
    number: ``For input string: "<s>"`` (observed)."""
    return JavaException("java.lang.NumberFormatException",
                         'For input string: "' + s + '"')


def utilities_exit(message: str) -> JavaTermination:
    """``blbutil.Utilities.exit(String)``: the message and a line separator
    on standard error, exit status 1."""
    return JavaTermination(stderr=message + "\n", status=1)


def utilities_exit_with_cause(message: str,
                              cause: JavaException) -> JavaTermination:
    """``blbutil.Utilities.exit(String, Throwable)``.

    Java writes the cause's stack trace, the cause, the message and
    ``terminating program.`` to standard error and exits with status 1.  The
    stack trace's first line (the cause) is reproduced; its ``at ...`` frame
    lines are not.
    """
    t = cause.to_string()
    return JavaTermination(
        stderr=t + "\n" + t + "\n" + message + "\n" + "terminating program.\n",
        status=1)


def utilities_exit_throwable(cause: JavaException) -> JavaTermination:
    """``blbutil.Utilities.exit(Throwable)``.

    Java writes ``ERROR: `` and the cause, the stack trace and
    ``terminating program.`` to standard output and exits with status 1.  The
    stack trace's first line is reproduced; its frame lines are not.
    """
    t = cause.to_string()
    return JavaTermination(
        stdout="ERROR: " + t + "\n" + t + "\n" + "terminating program.\n",
        status=1)


# ---------------------------------------------------------------------------
# Text: UTF-16 code units and output encoding
# ---------------------------------------------------------------------------

_NON_BMP = re.compile("[\U00010000-\U0010FFFF]")
_SURROGATE = re.compile("[\ud800-\udfff]")
_SURROGATE_PAIR = re.compile("[\ud800-\udbff][\udc00-\udfff]")


def _split_surrogates(match: "re.Match[str]") -> str:
    c = ord(match.group()) - 0x10000
    return chr(0xD800 + (c >> 10)) + chr(0xDC00 + (c & 0x3FF))


def _join_surrogates(match: "re.Match[str]") -> str:
    hi, lo = match.group()
    return chr(0x10000 + ((ord(hi) - 0xD800) << 10) + (ord(lo) - 0xDC00))


def to_java_chars(text: str) -> str:
    """Represents ``text`` as Java does: characters outside the Basic
    Multilingual Plane become two UTF-16 surrogate code units."""
    if _NON_BMP.search(text) is None:
        return text
    return _NON_BMP.sub(_split_surrogates, text)


def java_encode(text: str) -> bytes:
    """Encodes a Java string as UTF-8, with a surrogate pair becoming one
    character and an unpaired surrogate becoming ``?`` (the replacement of
    Java's UTF-8 encoder)."""
    if _SURROGATE.search(text) is not None:
        text = _SURROGATE_PAIR.sub(_join_surrogates, text)
        text = _SURROGATE.sub("?", text)
    return text.encode("utf-8")


def java_file_string(path: str) -> str:
    """How ``java.io.File`` shows a path given on the command line on Unix:
    repeated ``/`` are collapsed and a trailing ``/`` is removed."""
    normalized = re.sub("/+", "/", path)
    if len(normalized) > 1 and normalized.endswith("/"):
        normalized = normalized[:-1]
    return normalized


# ---------------------------------------------------------------------------
# Characters and strings (Java SE API documentation)
# ---------------------------------------------------------------------------

# String.trim(): "space is defined as any character whose codepoint is less
# than or equal to 'U+0020'".
_TRIM_CHARS = "".join(chr(code) for code in range(0x21))

# Character.isWhitespace(char): the Javadoc's list of control characters,
# and the Unicode space, line and paragraph separators except the three
# non-breaking spaces.
_WHITESPACE_CONTROLS = frozenset("\t\n\x0b\x0c\r\x1c\x1d\x1e\x1f")
_NON_BREAKING = frozenset("   ")


def java_trim(s: str) -> str:
    """``String.trim()``."""
    return s.strip(_TRIM_CHARS)


def java_is_whitespace(c: str) -> bool:
    """``Character.isWhitespace(char)``."""
    if c in _WHITESPACE_CONTROLS:
        return True
    return (c not in _NON_BREAKING
            and unicodedata.category(c) in ("Zs", "Zl", "Zp"))


def java_is_digit(c: str) -> bool:
    """``Character.isDigit(char)``: Unicode general category Nd."""
    return unicodedata.category(c) == "Nd"


def java_digit(c: str) -> int:
    """``Character.digit(char, 10)``: the value of a decimal digit character
    (category Nd), otherwise -1.  A surrogate code unit is not a digit."""
    return unicodedata.decimal(c, -1)


def java_is_acgtn(c: str) -> bool:
    """``Character.toUpperCase(c)`` is one of A, C, G, T, N.

    No character outside ASCII has one of these letters as its upper-case
    form (checked against every Basic Multilingual Plane character with
    Python's Unicode 15.1 tables), so the test reduces to the ASCII letters.
    """
    return c in "ACGTNacgtn"


def java_substring(s: str, begin: int, end: Optional[int] = None) -> str:
    """``String.substring(begin, end)``, which the jar calls to quote the
    first 80 characters of a VCF record.  Out-of-range bounds raise
    StringIndexOutOfBoundsException with the message Java 25 gives
    (observed): ``Range [begin, end) out of bounds for length n``."""
    if end is None:
        end = len(s)
    if begin < 0 or begin > end or end > len(s):
        raise JavaException(
            "java.lang.StringIndexOutOfBoundsException",
            "Range [" + str(begin) + ", " + str(end)
            + ") out of bounds for length " + str(len(s)))
    return s[begin:end]


def java_array_to_string(items: List[str]) -> str:
    """``java.util.Arrays.toString(Object[])``: ``[a, b, c]``."""
    return "[" + ", ".join(items) + "]"


# ---------------------------------------------------------------------------
# blbutil.StringUtil (translated from Browning's Java)
# ---------------------------------------------------------------------------

_JAVA_TOKEN = re.compile("[^\x00-\x20]+")


def get_fields_ws(s: str) -> List[str]:
    """``StringUtil.getFields(String)``: the fields of the trimmed string,
    separated by runs of characters whose code is at most U+0020."""
    return _JAVA_TOKEN.findall(s)


def get_fields_ws_limit(s: str, limit: int) -> List[str]:
    """``StringUtil.getFields(String, int)``: as ``get_fields_ws`` but with
    at most ``limit`` fields; the last field is the rest of the trimmed
    string."""
    if limit < 2:
        raise illegal_argument("limit: " + str(limit))
    tokens = list(_JAVA_TOKEN.finditer(s))
    if len(tokens) <= limit:
        return [m.group() for m in tokens]
    fields = [m.group() for m in tokens[:limit - 1]]
    fields.append(java_trim(s[tokens[limit - 1].start():]))
    return fields


def get_fields_delim(s: str, delimiter: str) -> List[str]:
    """``StringUtil.getFields(String, char)``: every field between
    delimiters, empty fields included."""
    return s.split(delimiter)


def get_fields_delim_limit(s: str, delimiter: str, limit: int) -> List[str]:
    """``StringUtil.getFields(String, char, int)``: at most ``limit``
    fields; the last field is the rest of the string."""
    if limit < 2:
        raise illegal_argument("limit: " + str(limit))
    return s.split(delimiter, limit - 1)


# ---------------------------------------------------------------------------
# Integer.parseInt (Java SE API documentation)
# ---------------------------------------------------------------------------

_ASCII_INT = re.compile("[+-]?[0-9]{1,9}\\Z")


def java_parse_int(s: str) -> int:
    """``Integer.parseInt(String)``.

    The Javadoc: every character must be a decimal digit
    (``Character.digit(c, 10)`` is not -1), except that the first may be
    ``-`` or ``+`` when the string is longer than one character, and the
    value must fit in an int.  Otherwise NumberFormatException, with the
    message ``For input string: "<s>"`` (observed).
    """
    if _ASCII_INT.match(s):
        return int(s)
    negative = s[:1] == "-"
    body = s[1:] if s[:1] in ("-", "+") else s
    values = [java_digit(c) for c in body]
    if not values or min(values) < 0:
        raise number_format_for_input(s)
    while len(values) > 1 and values[0] == 0:
        del values[0]
    if len(values) > 10:
        # More than 10 significant digits is outside the int range; this
        # also keeps Python's limit on long integer strings out of reach.
        raise number_format_for_input(s)
    magnitude = int("".join(str(v) for v in values))
    value = -magnitude if negative else magnitude
    if not INT_MIN <= value <= INT_MAX:
        raise number_format_for_input(s)
    return value


def to_int32(value: int) -> int:
    """Wraps an integer to 32 bits, as Java ``int`` arithmetic does."""
    return ((value + (1 << 31)) & 0xFFFFFFFF) - (1 << 31)


# ---------------------------------------------------------------------------
# Double.parseDouble and Float.parseFloat (Java SE API documentation)
# ---------------------------------------------------------------------------
#
# The Javadoc of Double.valueOf(String) defines the accepted strings: after
# leading and trailing characters up to U+0020 are removed (as by
# String.trim), an optional sign followed by "NaN", "Infinity", or a decimal
# or hexadecimal floating-point literal of the Java Language Specification
# (section 3.10.2) written without underscores, which may end with one of
# the type suffixes f, F, d, D.  The result is the exact value of the literal
# rounded to the nearest double (Double.parseDouble) or float
# (Float.parseFloat), ties to even.
#
# The two messages of the NumberFormatException are observed from the jar
# (tests/parity, float_string cases): "multiple points" when the run of
# digits and points at the start of the unsigned string (of hexadecimal
# digits and points after a 0x or 0X prefix) contains two points, and
# 'For input string: "<string>"' otherwise, with a long string shortened as
# described at _MESSAGE_LIMIT.  The empty string is never parsed by this
# program, so its message is not observed.

_FLOAT_LITERAL = re.compile(
    "(?P<sign>[+-]?)(?:"
    "(?P<nan>NaN)"
    "|(?P<inf>Infinity)"
    "|(?P<dec>(?:[0-9]+(?:\\.[0-9]*)?|\\.[0-9]+)(?:[eE][+-]?[0-9]+)?)[fFdD]?"
    "|0[xX](?P<hexint>[0-9a-fA-F]*)(?:\\.(?P<hexfrac>[0-9a-fA-F]*))?"
    "[pP](?P<binexp>[+-]?[0-9]+)[fFdD]?"
    ")\\Z")
_TWO_POINTS = re.compile(
    "[+-]?(?:0[xX][0-9a-fA-F]*\\.[0-9a-fA-F]*\\.|[0-9]*\\.[0-9]*\\.)")
# The jar shortens a long string in its NumberFormatException message to
# its first 497 and last 498 characters around " ... " (observed for a
# 1501-character string; tests/parity float_string cases record the
# lengths around the threshold).
_MESSAGE_LIMIT = 1000
_MESSAGE_HEAD = 497
_MESSAGE_TAIL = 498
_PLAIN_DECIMAL = re.compile("-?[0-9]{1,15}(?:\\.[0-9]{1,15})?\\Z")

# Beyond these binary magnitudes a value is infinite or zero in both float
# and double; the limits keep exact arithmetic small.
_MAX_BINARY_MAGNITUDE = 1100
_MIN_BINARY_MAGNITUDE = -1200


class _FloatLiteral:
    """An accepted floating-point string: its sign and either a special
    value (``nan``, ``inf``), the decimal text for Python's ``float``, or the
    exact value of a hexadecimal literal (``None`` meaning out of range)."""

    __slots__ = ("negative", "special", "decimal", "exact")

    def __init__(self, negative: bool, special: str = "",
                 decimal: str = "", exact: Optional[Fraction] = None):
        self.negative = negative
        self.special = special
        self.decimal = decimal
        self.exact = exact


def _read_float_literal(s: str) -> _FloatLiteral:
    text = java_trim(s)
    m = _FLOAT_LITERAL.match(text)
    if m is None or (m.group("binexp") is not None
                     and not (m.group("hexint") or m.group("hexfrac"))):
        if _TWO_POINTS.match(text):
            raise JavaException("java.lang.NumberFormatException",
                                "multiple points")
        if len(s) > _MESSAGE_LIMIT:
            s = s[:_MESSAGE_HEAD] + " ... " + s[len(s) - _MESSAGE_TAIL:]
        raise number_format_for_input(s)
    negative = m.group("sign") == "-"
    if m.group("nan"):
        return _FloatLiteral(negative, special="nan")
    if m.group("inf"):
        return _FloatLiteral(negative, special="inf")
    if m.group("dec") is not None:
        return _FloatLiteral(negative, decimal=m.group("dec"))
    digits = m.group("hexint") + (m.group("hexfrac") or "")
    significand = int(digits, 16)
    if significand == 0:
        return _FloatLiteral(negative, exact=Fraction(0))
    binexp = m.group("binexp")
    exp_digits = binexp.lstrip("+-").lstrip("0") or "0"
    if len(exp_digits) > 12:
        # |exponent| of at least 10**12: infinite or zero whatever the
        # significand's length (which a field could not approach).
        return _FloatLiteral(negative, special="" if binexp[0] == "-"
                             else "inf", exact=Fraction(0))
    exponent = int(binexp) - 4 * len(m.group("hexfrac") or "")
    magnitude = exponent + significand.bit_length()
    if magnitude > _MAX_BINARY_MAGNITUDE:
        return _FloatLiteral(negative, special="inf")
    if magnitude < _MIN_BINARY_MAGNITUDE:
        return _FloatLiteral(negative, exact=Fraction(0))
    return _FloatLiteral(negative, exact=Fraction(significand)
                         * Fraction(2) ** exponent)


def java_parse_double(s: str) -> float:
    """``Double.parseDouble(String)``, correctly rounded (ties to even)."""
    if _PLAIN_DECIMAL.match(s):
        return float(s)
    lit = _read_float_literal(s)
    sign = -1.0 if lit.negative else 1.0
    if lit.special == "nan":
        return math.nan
    if lit.special == "inf":
        return sign * math.inf
    if lit.decimal:
        return sign * float(lit.decimal)
    try:
        return sign * float(lit.exact)
    except OverflowError:
        return sign * math.inf


def _round_to_float32(q: Fraction) -> float:
    """Rounds a non-negative exact value to the nearest float32 value (ties
    to even; values from halfway above the largest float round to
    infinity), returned as a Python float."""
    if q == 0:
        return 0.0
    e = q.numerator.bit_length() - q.denominator.bit_length()
    if Fraction(2) ** e > q:
        e -= 1
    if e > 127:
        return math.inf
    quantum = Fraction(2) ** max(e - 23, -149)
    n = round(q / quantum)          # Fraction rounding is half-to-even
    value = n * quantum
    if value >= 2 ** 128:
        return math.inf
    return float(value)


def _is_float32_halfway(d: float) -> bool:
    """True if the positive double ``d`` lies exactly halfway between two
    adjacent float32 values.  Only there can rounding the double to float32
    differ from rounding the original decimal value."""
    m, e = math.frexp(d)
    significand = int(m * 9007199254740992.0)   # m * 2**53, exact
    dropped = max((e - 1) - 23, -149) - (e - 53)
    if dropped <= 0 or dropped > 60:
        return False
    return significand % (1 << dropped) == 1 << (dropped - 1)


def _double_to_float32(d: float) -> float:
    """Rounds a finite double to float32 (ties to even)."""
    try:
        return struct.unpack("<f", struct.pack("<f", d))[0]
    except OverflowError:
        return math.copysign(math.inf, d)


def _decimal_fraction(text: str) -> Fraction:
    """The exact value of an unsigned decimal literal such as ``1.25e-3``
    (through ``decimal.Decimal``, which has no limit on the number of
    digits)."""
    return Fraction(Decimal(text))


def java_parse_float(s: str) -> float:
    """``Float.parseFloat(String)``: rounded to float32 directly from the
    exact value of the string (not through double), returned as a Python
    float."""
    if _PLAIN_DECIMAL.match(s):
        negative = s.startswith("-")
        text = s.lstrip("-")
    else:
        lit = _read_float_literal(s)
        negative = lit.negative
        if lit.special == "nan":
            return math.nan
        if lit.special == "inf":
            return -math.inf if negative else math.inf
        if not lit.decimal:
            value = _round_to_float32(lit.exact)
            return -value if negative else value
        text = lit.decimal
    d = float(text)
    if d == 0.0 or math.isinf(d):
        value = d
    elif _is_float32_halfway(d):
        value = _round_to_float32(_decimal_fraction(text))
    else:
        value = _double_to_float32(d)
    return -value if negative else value


# ---------------------------------------------------------------------------
# Ordering: Float.compare / Double.compare and Arrays.binarySearch
# ---------------------------------------------------------------------------

def java_compare_key(x: float) -> Tuple[int, float, float]:
    """A sort key for the total order of ``Float.compare`` and
    ``Double.compare`` (Javadoc): -0.0 is less than 0.0, and NaN is greater
    than every other value and equal to itself."""
    if x != x:
        return (1, 0.0, 0.0)
    return (0, x, math.copysign(1.0, x))


def _bisect(n: int, compare) -> int:
    """Binary search over indices 0..n-1 by repeated halving: each probe is
    the lower middle of the remaining inclusive range.  ``compare(i)`` is
    negative, zero or positive as element i is below, at or above the key.
    Returns the index of a match, or -(insertion point) - 1."""
    low, high = 0, n - 1
    while low <= high:
        probe = (low + high) // 2
        c = compare(probe)
        if c == 0:
            return probe
        if c < 0:
            low = probe + 1
        else:
            high = probe - 1
    return -(low + 1)


def java_binary_search_int(a: List[int], key: int) -> int:
    """``Arrays.binarySearch(int[], int)``."""
    return _bisect(len(a), lambda i: (a[i] > key) - (a[i] < key))


def java_binary_search_double(a: List[float], key: float) -> int:
    """``Arrays.binarySearch(double[], double)``, ordered as by
    ``Double.compare``.  The Javadoc leaves open which of several equal
    elements is found; this one returns the one the halving reaches first,
    which agrees with the jar where the genetic map repeats the cM value
    that the 5 cM end rule searches for (tests/parity, map_equal_cm_start,
    map_equal_cm_end, map_equal_cm_many)."""
    k = java_compare_key(key)

    def compare(i: int) -> int:
        v = java_compare_key(a[i])
        return (v > k) - (v < k)

    return _bisect(len(a), compare)


# ---------------------------------------------------------------------------
# DecimalFormat("0.##") and DecimalFormat("0.###")
# ---------------------------------------------------------------------------

def java_decimal_format(x: float, maximum_fraction_digits: int) -> str:
    """``new DecimalFormat("0.##").format(double)`` (2 fraction digits) or
    ``"0.###"`` (3), as the jar prints.

    From the DecimalFormat and RoundingMode.HALF_EVEN documentation, checked
    against the jar's output for about 97,000 values (tests/parity,
    vectors_decimal_format; Java 25.0.4.1):

    - NaN prints ``NaN``; infinities print ``∞`` with the sign.
    - A value is rounded to the pattern's fraction digits, ties to even,
      using its exact binary value.  At least one integer digit is printed;
      there is no grouping and no trailing fraction zero.
    - A negative value, including -0.0, keeps its ``-`` even when it prints
      as zero (``-0``).
    - A value whose shortest decimal form has no more fraction digits than
      the pattern allows prints that form (Python's ``repr``) padded with
      zeros, not its rounded exact binary value; the jar does the same,
      including for values from 2**43 to 2**53 (observed).
    - For some integers of 2**53 and above the jar prints other digits (for
      example 2**60 as ``1152921504606846980`` and 1e23 as
      ``99999999999999990000000``); both forms convert back to the same
      number.  This conversion does not reproduce the jar's digits there: in
      the exact-conversion test 169 of 1,794 such values printed
      differently, the smallest about 1.84e16.
    """
    if x != x:
        return "NaN"
    sign = "-" if math.copysign(1.0, x) < 0 else ""
    if math.isinf(x):
        return sign + "∞"
    a = abs(x)
    shortest = Decimal(repr(a))
    if shortest.as_tuple().exponent >= -maximum_fraction_digits:
        value = shortest
    else:
        with localcontext() as context:
            context.prec = 1000
            value = Decimal(a).quantize(
                Decimal(1).scaleb(-maximum_fraction_digits),
                rounding=ROUND_HALF_EVEN)
    text = format(value, "f")
    if "." in text:
        text = text.rstrip("0").rstrip(".")
    return sign + text
