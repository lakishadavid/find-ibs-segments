# find-ibs-segments.  This module comes from merge-ibs-segments, a modified
# Python translation of merge-ibd-segments from Refined-IBD version 17Jan20.102.
#
# Copyright (C) 2014-2016 Brian L. Browning
# Copyright (c) 2009 The Broad Institute
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
# Python from src/blbutil/InputIt.java (including where it ends the program
# through Utilities.exit, which is translated in java_compat.py) in
# refined-ibd.17Jan20.102.zip (SHA-256
# 0b0abf48528ec53d3fec7f21a9742c6e5960857c5a225b0e49346179bc45e6ea). Later
# changes are recorded in the Git history.
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
#
# _is_valid_bgzf_header, _bgzf_chunks, their constants and their error
# messages are translated from src/net/sf/samtools/util/
# BlockCompressedInputStream.java, BlockCompressedStreamConstants.java and
# BlockGunzipper.java (Picard's sam-jdk, bundled in the same archive).
# Those files carry this notice, reproduced unchanged:
#
# The MIT License
#
# Copyright (c) 2009 The Broad Institute
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
"""Line input, converted from blbutil.InputIt.

``InputIt.fromGzipFile(File)`` opens a file, decompresses it when its name
ends in ``.gz`` (BGZF through net.sf.samtools.util.BlockCompressedInputStream,
other gzip files through java.util.zip.GZIPInputStream), and returns its
lines as java.io.BufferedReader.readLine() does: a line ends at ``\\n``,
``\\r`` or ``\\r\\n`` (Javadoc).  ``InputIt.fromStdIn()`` does the same for
standard input without decompression.  Bytes are decoded as UTF-8, the
default charset of Java 18 and later, with malformed input replaced by
U+FFFD.  Lines are produced lazily, so an error in the middle of a file is
reported when that point is reached, as in Java.

A read error ends the program through ``Utilities.exit("Error reading " +
reader, e)``.  InputIt reads the first line in its constructor, through a
BufferedReader of 4,194,304 characters, and names the decompressing stream
in that message; later reads name the BufferedReader (observed).  Java also
prints the object's identity hash (``@76ed5528``), which is not reproduced.
This conversion names the stream for an error in the first 4,194,304
characters, or before the end of the first line if that is longer; for BGZF
files, where Java's first read can stop at a block boundary, this boundary
has not been checked against the jar, and neither have the messages for a
damaged BGZF file.
"""

from __future__ import annotations

import codecs
import os
import re
import struct
import zlib
from typing import BinaryIO, Iterator, List

from .java_compat import (JavaException, JavaTermination, java_file_string,
                          to_java_chars, utilities_exit_with_cause)

_CHUNK = 1 << 20
_FIRST_READ_CHARS = 1 << 22
_BGZF_HEADER_LENGTH = 18
_LINE_TERMINATOR = re.compile("\r\n|\r|\n")

GZIP_STREAM = "java.util.zip.GZIPInputStream"
BGZF_STREAM = "net.sf.samtools.util.BlockCompressedInputStream"
FILE_STREAM = "java.io.FileInputStream"
READER = "java.io.BufferedReader"


class _ReadFailure(Exception):
    """A decompression error; the reader named in the message depends on
    how far reading had come."""

    def __init__(self, cause: JavaException):
        super().__init__(cause.to_string())
        self.cause = cause


def _is_valid_bgzf_header(header: bytes) -> bool:
    """``BlockCompressedInputStream.isValidBlockHeader``: gzip magic bytes,
    the FEXTRA flag, an extra field of length 6 and the ``BC`` subfield."""
    return (len(header) == _BGZF_HEADER_LENGTH
            and header[0] == 31 and header[1] == 139
            and (header[3] & 4) != 0
            and header[10] == 6
            and header[12] == 66 and header[13] == 67)


def _zlib_cause(error: Exception) -> JavaException:
    """The Java exception for a decompression failure.  zlib's own message
    follows the last colon of Python's message; Java's Inflater reports the
    same zlib message in a ZipException."""
    text = str(error)
    if "incorrect data check" in text or "incorrect length check" in text:
        return JavaException("java.util.zip.ZipException",
                             "Corrupt GZIP trailer")
    return JavaException("java.util.zip.ZipException",
                         text.rsplit(": ", 1)[-1])


def _plain_chunks(f: BinaryIO) -> Iterator[bytes]:
    while True:
        chunk = f.read(_CHUNK)
        if not chunk:
            return
        yield chunk


def _gzip_chunks(f: BinaryIO) -> Iterator[bytes]:
    """The members of a gzip file, one after another.  Bytes after a member
    that do not begin another gzip member are ignored."""
    buffered = b""
    eof = False

    def fill(minimum: int) -> None:
        nonlocal buffered, eof
        while len(buffered) < minimum and not eof:
            chunk = f.read(_CHUNK)
            if not chunk:
                eof = True
            buffered += chunk

    first = True
    while True:
        fill(3)
        if not first and (len(buffered) < 3 or buffered[0] != 31
                          or buffered[1] != 139 or buffered[2] != 8):
            return
        first = False
        d = zlib.decompressobj(31)
        data = buffered
        buffered = b""
        while True:
            try:
                out = d.decompress(data)
            except zlib.error as error:
                raise _ReadFailure(_zlib_cause(error))
            if out:
                yield out
            if d.eof:
                buffered = d.unused_data
                break
            data = f.read(_CHUNK)
            if not data:
                raise _ReadFailure(JavaException(
                    "java.io.EOFException",
                    "Unexpected end of ZLIB input stream"))
        fill(1)
        if not buffered:
            return


def _bgzf_chunks(f: BinaryIO) -> Iterator[bytes]:
    """BGZF blocks.  Block CRCs are not checked, as in
    BlockCompressedInputStream's default setting."""
    while True:
        header = f.read(_BGZF_HEADER_LENGTH)
        if not header:
            return
        if len(header) != _BGZF_HEADER_LENGTH:
            raise _ReadFailure(JavaException("java.io.IOException",
                                             "Premature end of file"))
        if not _is_valid_bgzf_header(header):
            raise _ReadFailure(JavaException("java.io.IOException",
                                             "Invalid GZIP header"))
        block_size = struct.unpack("<H", header[16:18])[0] + 1
        rest = f.read(block_size - _BGZF_HEADER_LENGTH)
        if len(rest) != block_size - _BGZF_HEADER_LENGTH:
            raise _ReadFailure(JavaException("java.io.IOException",
                                             "Premature end of file"))
        try:
            out = zlib.decompress(rest[:-8], -15)
        except zlib.error as error:
            raise _ReadFailure(_zlib_cause(error))
        if out:
            yield out


def _split_lines(text: str) -> List[str]:
    """Splits at ``\\r\\n``, ``\\r`` and ``\\n``; the last element is the
    text after the last terminator (possibly empty)."""
    return _LINE_TERMINATOR.split(text)


def java_lines_from_chunks(chunks: Iterator[bytes],
                           stream_name: str) -> Iterator[str]:
    """The lines ``BufferedReader.readLine()`` returns, from byte chunks.
    ``stream_name`` is the Java class named in an early read error."""
    decoder = codecs.getincrementaldecoder("utf-8")(errors="replace")
    pending = ""
    decoded = 0
    first_line_done = False

    def failure(error: _ReadFailure) -> JavaTermination:
        early = decoded < _FIRST_READ_CHARS or not first_line_done
        return utilities_exit_with_cause(
            "Error reading " + (stream_name if early else READER),
            error.cause)

    try:
        for chunk in chunks:
            new_text = decoder.decode(chunk)
            decoded += len(new_text)
            text = pending + new_text
            keep = ""
            if text.endswith("\r"):
                # A following "\n" belongs to the same line terminator.
                text = text[:-1]
                keep = "\r"
            parts = _split_lines(text)
            pending = parts.pop() + keep
            for part in parts:
                first_line_done = True
                yield to_java_chars(part)
    except _ReadFailure as error:
        raise failure(error)
    text = pending + decoder.decode(b"", final=True)
    if text:
        parts = _split_lines(text)
        if parts[-1] == "":
            parts.pop()
        for part in parts:
            yield to_java_chars(part)


def from_gzip_file(path: str) -> Iterator[str]:
    """``InputIt.fromGzipFile(new File(path))``.

    Raises JavaTermination at once if the file cannot be opened or, for a
    ``.gz`` name that is not BGZF, if the gzip header is not valid; later
    errors are raised while lines are read.
    """
    display = java_file_string(path)
    try:
        f = open(path, "rb")
    except OSError as error:
        reason = error.strerror or str(error)
        cause = JavaException("java.io.FileNotFoundException",
                              display + " (" + reason + ")")
        raise utilities_exit_with_cause("Error opening " + display, cause)
    if not os.path.basename(path).endswith(".gz"):
        return _generate(f, _plain_chunks(f), FILE_STREAM)
    header = f.read(_BGZF_HEADER_LENGTH)
    f.seek(0)
    if _is_valid_bgzf_header(header):
        return _generate(f, _bgzf_chunks(f), BGZF_STREAM)
    cause = None
    if len(header) < 2:
        cause = JavaException("java.io.EOFException")
    elif header[0] != 31 or header[1] != 139:
        cause = JavaException("java.util.zip.ZipException",
                              "Not in GZIP format")
    elif len(header) < 3:
        cause = JavaException("java.io.EOFException")
    elif header[2] != 8:
        cause = JavaException("java.util.zip.ZipException",
                              "Unsupported compression method")
    if cause is not None:
        f.close()
        raise utilities_exit_with_cause("Error reading " + display, cause)
    return _generate(f, _gzip_chunks(f), GZIP_STREAM)


def _generate(f: BinaryIO, chunks: Iterator[bytes],
              stream_name: str) -> Iterator[str]:
    try:
        yield from java_lines_from_chunks(chunks, stream_name)
    finally:
        f.close()


def from_std_in(stream: BinaryIO) -> Iterator[str]:
    """``InputIt.fromStdIn()``: the lines of standard input."""
    return java_lines_from_chunks(_plain_chunks(stream), "java.io.InputStream")
