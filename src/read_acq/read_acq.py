"""Functions and classes for reading and writing the .acq format."""

from __future__ import annotations

import contextlib
import io
import re
import warnings
from collections.abc import Iterable, Iterator
from pathlib import Path
from typing import BinaryIO, ClassVar

import attrs
import numpy as np
import tqdm

from .codec import _decode_into, _decode_line, _encode


class ACQError(Exception):
    """Base class for errors in the ACQ module."""


class ACQLineError(ACQError):
    """Error in parsing a line of an ACQ file."""


@attrs.define()
class CommentLine:
    """A class for storing a comment line."""

    swpos: int = attrs.field(converter=int)
    adcmax: float = attrs.field(converter=float)
    adcmin: float = attrs.field(converter=float)
    temp: float = attrs.field(converter=float)
    nblk: int = attrs.field(converter=int)
    nspec: int = attrs.field(converter=int)
    resolution: float = attrs.field(
        converter=attrs.converters.optional(float), default=None
    )
    data_drops: int = attrs.field(
        converter=attrs.converters.optional(int), default=None
    )

    _sw = r"swpos (?P<swpos>\d)"
    _res = r"resolution \s*(?P<resolution>\d+(\.\d*)?|\.\d+)"
    _acmax = r"adcmax \s*(?P<adcmax>[-+]?\d+(\.\d*)?|\.\d+)"
    _acmin = r"adcmin \s*(?P<adcmin>[-+]?\d+(\.\d*)?|\.\d+)"
    _temp = r"temp \s*(?P<temp>\d+) C"
    _nblk = r"nblk \s*(?P<nblk>\d+)"
    _nspec = r"nspec \s*(?P<nspec>\d+)"
    _dd = r"data_drops \s*(?P<data_drops>\d+)"

    pxspec = re.compile(f"# {_sw} {_res} {_acmax} {_acmin} {_temp} {_nblk} {_nspec}")
    fastspec = re.compile(f"# {_sw} {_dd} {_acmax} {_acmin} {_temp} {_nblk} {_nspec}")

    @classmethod
    def read(cls, line, fastspec: bool | None = None):
        """Read an ACQ comment line as a CommentLine object."""
        line = line.strip()
        if fastspec is False:
            match = re.match(cls.pxspec, line)
        elif fastspec is True:
            match = re.match(cls.fastspec, line)
        else:
            match = re.match(cls.pxspec, line)
            if match is None:
                match = re.match(cls.fastspec, line)
        if match is None:
            raise ACQError(f"Could not parse line: '{line}'")

        return cls(**match.groupdict())


@attrs.define()
class DataLine:
    """A class for storing a data line."""

    time: str = attrs.field()
    swpos: int = attrs.field(converter=int)
    freqmin: float = attrs.field(converter=float)
    deltaf: float = attrs.field(converter=float)
    freqmax: float = attrs.field(converter=float)
    thing: float = attrs.field(converter=float)
    spectrum: np.ndarray | None = attrs.field(default=None)

    _time = r"(?P<time>\d{4}:\d{3}:\d{2}:\d{2}:\d{2})"
    _swpos = r"(?P<swpos>\d{1})"
    _float = r"[-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?"
    _fmin = f"(?P<freqmin>{_float})"
    _df = f"(?P<deltaf>{_float})"
    _fmax = f"(?P<freqmax>{_float})"
    _last = f"(?P<thing>{_float})"

    regex = re.compile(rf"{_time} {_swpos} \s*{_fmin} \s*{_df} \s*{_fmax} \s*{_last}")

    @classmethod
    def read(cls, line, read_spectrum: bool = True):
        """Read an ACQ data line as a DataLine object."""
        try:
            front, back = line.split(" spectrum ")
        except ValueError:
            raise ACQLineError(
                f"Could not parse line: '{line[:100]}' -- probably incomplete"
            ) from None

        match = re.match(cls.regex, front)
        if match is None:
            raise ACQError(f"Could not parse line front-matter: {front}")

        spec = _decode_line(back.lstrip()) if read_spectrum else None
        return cls(spectrum=spec, **match.groupdict())


@attrs.define
class DataEntry:
    """A class for storing a data entry."""

    comment: CommentLine = attrs.field()
    data: DataLine = attrs.field()

    @data.validator
    def _check_data(self, attribute, value):
        if value.swpos != self.comment.swpos:
            raise ACQLineError("swpos of comment and data do not match")
        if value.spectrum is not None and value.spectrum.shape[0] != self.comment.nspec:
            raise ACQLineError("nspec and length of spectrum do not match")

    @classmethod
    def read(cls, lines, read_spectrum: bool = True) -> DataEntry:
        """Read an ACQ data entry (two lines) as a DataEntry object."""
        comment = CommentLine.read(lines[0])
        data = DataLine.read(lines[1], read_spectrum=read_spectrum)
        return cls(comment, data)


class Ancillary:
    """The ancillary data of an ACQ file."""

    header_char = ";"

    _splits = re.compile(r"[\d\.:]+")
    DTYPES: ClassVar = {
        "adcmax": np.float32,
        "adcmin": np.float32,
        "times": "S17",
    }

    def __init__(self, fname: str | Path):
        fname = Path(fname)
        with fname.open("rb", buffering=0) as fl:
            self._init(fname, _read_head(_Reader(fl)))

    @classmethod
    def _from_reader(cls, fname: Path, reader: _Reader) -> Ancillary:
        """Create the ancillary data from a file that is already open.

        The reader is left at the start of the file, keeping the bytes it has read.
        """
        anc = cls.__new__(cls)
        anc._init(fname, _read_head(reader))
        return anc

    def _init(self, fname: Path, lines: list[str]):
        self.fastspec_version = self._fastspec_version_from(lines)
        self.meta = self._metadata_from(fname, lines)

        self.data = {
            "adcmax": [],  # np.zeros((self.size, 3), dtype=np.float32),
            "adcmin": [],  # np.zeros((self.size, 3), dtype=np.float32),
            "times": [],  # np.zeros((self.size, 3), dtype="S17"),
        }
        if "data_drops" in self.meta:
            self.data["data_drops"] = []  # np.zeros((self.size, 3), dtype=int)

    @staticmethod
    def _head_lines(fname: Path) -> list[str]:
        with fname.open("rb", buffering=0) as fl:
            return _read_head(_Reader(fl))

    def _get_fastspec_version(self, fname: Path):
        return self._fastspec_version_from(self._head_lines(fname))

    def _fastspec_version_from(self, lines: list[str]):
        first_line = lines[0] if lines else ""
        if first_line.startswith(self.header_char):
            return first_line.split("FASTSPEC")[-1]
        else:
            return None

    def _read_header(self, fname: Path):
        return self._header_from(fname, self._head_lines(fname))

    def _header_from(self, fname: Path, lines: list[str]):
        out = {}

        name_pattern = re.compile(r"[a-zA-Z_]+")
        type_order = [int, float, str]

        for line in lines:
            if not line.startswith(self.header_char):
                break

            if line.startswith("; FASTSPEC"):
                name = "fastspec_version"
                val = line.split()[-1]
            else:
                try:
                    name, val = line.split(": ")
                except ValueError:
                    warnings.warn(
                        f"In file {fname}, item {line} has no value", stacklevel=1
                    )
                    name = line.split(":")[0]
                    val = ""

                name = name_pattern.findall(name)[0]

            for tp in type_order:
                try:
                    out[name] = tp(val.split()[0])
                    break
                except IndexError:
                    with contextlib.suppress(ValueError):
                        out[name] = tp(val)
                except ValueError:
                    pass

        return out

    def read_metadata(self, fname: Path):
        """Read the metadata of the ACQ file."""
        return self._metadata_from(fname, self._head_lines(fname))

    def _metadata_from(self, fname: Path, lines: list[str]):
        out = self._header_from(fname, lines)

        it = iter(lines)
        for line in it:
            if line.startswith("#"):
                comment = CommentLine.read(line)
                data = DataLine.read(next(it), read_spectrum=False)
                break
        else:
            raise ACQError(f"No comment line found in file {fname}.")

        out.update(
            {
                "temperature": comment.temp,
                "nblk": comment.nblk,
                "nfreq": comment.nspec,
                "freq_min": data.freqmin,
                "freq_max": data.freqmax,
                "freq_res": data.deltaf,
            }
        )
        if comment.resolution is not None:
            out["resolution"] = comment.resolution
        if comment.data_drops is not None:
            out["data_drops"] = comment.data_drops

        return out

    @property
    def frequencies(self):
        """The frequencies associated with the spectrum measurements."""
        df = self.meta["freq_max"] / self.meta["nfreq"]

        # See edges-cal.tools.EdgesFrequencyRange for justification of using this
        # form.
        return np.arange(self.meta["freq_min"], self.meta["freq_max"], df)

    def append(self, datas: tuple[DataEntry, DataEntry, DataEntry]):
        """Append a set of data to the ancillary data."""
        if tuple(d.data.swpos for d in datas) != (0, 1, 2):
            raise ValueError(
                "swpos of datas must be (0, 1, 2), got "
                f"{tuple(d.data.swpos for d in datas)}"
            )

        self.data["adcmax"].append([d.comment.adcmax for d in datas])
        self.data["adcmin"].append([d.comment.adcmin for d in datas])
        self.data["times"].append([d.data.time for d in datas])

        if "data_drops" in self.meta:
            self.data["data_drops"].append([d.comment.data_drops for d in datas])

    def complete(self):
        """Convert the ancillary data to numpy arrays."""
        for key in self.data:
            self.data[key] = np.array(
                self.data[key], dtype=self.DTYPES.get(key, np.float32)
            )


_SPECTRUM_SEP = " spectrum "


def _encoded_length(line: str) -> int:
    """Return the number of encoded values in a data line, without decoding them."""
    start = line.index(_SPECTRUM_SEP) + len(_SPECTRUM_SEP)
    return len(line[start:].lstrip()) // 4


def _warn_nul_padded(fname: Path):
    warnings.warn(
        f"File {fname} contains NUL bytes; it was probably not fully "
        "written. Returning the complete cycles read so far.",
        stacklevel=2,
    )


# Bytes first read at the start of a file, for its header and first entry. This is
# doubled until they have all been read.
_HEAD_READ_SIZE = 4096

# Bytes read at a time when reading a whole line.
_LINE_READ_SIZE = 65536


class _Reader:
    """Read a binary file at a moving position, keeping the bytes read ahead of it.

    This gives control over how much is read, and where: each read of the file is
    a seek and a read of a given size.
    """

    def __init__(self, fl: BinaryIO):
        self.fl = fl
        self.pos = 0
        self.ahead = b""  # The bytes of the file from pos.
        self.eof = False  # Whether ``ahead`` reaches the end of the file.

    def read(self, pos: int, size: int) -> tuple[bytes, bool]:
        """Read ``size`` bytes at ``pos``, and whether that reaches the end of the file.

        This does not move the reader.
        """
        self.fl.seek(pos)
        data = self.fl.read(size)
        # An unbuffered read may return fewer bytes than asked for.
        while 0 < len(data) < size and (more := self.fl.read(size - len(data))):
            data += more
        return data, len(data) < size

    def move(self, pos: int, ahead: bytes, eof: bool):
        """Move to ``pos``, where ``ahead`` has already been read."""
        self.pos, self.ahead, self.eof = pos, ahead, eof

    def skip(self, n: int):
        """Move forward ``n`` bytes."""
        self.pos += n
        self.ahead = self.ahead[n:]

    def peek(self, size: int) -> bytes:
        """Return at least the next ``size`` bytes (fewer at the end of the file)."""
        if len(self.ahead) < size and not self.eof:
            more, self.eof = self.read(
                self.pos + len(self.ahead), size - len(self.ahead)
            )
            self.ahead += more
        return self.ahead

    def readline(self) -> bytes:
        """Return the line at the current position, without moving past it."""
        end = self.ahead.find(b"\n") + 1
        parts = [self.ahead]
        nread = len(self.ahead)
        while not end and not self.eof:
            more, self.eof = self.read(self.pos + nread, _LINE_READ_SIZE)
            if nl := more.find(b"\n") + 1:
                end = nread + nl
            parts.append(more)
            nread += len(more)
        self.ahead = b"".join(parts)
        return self.ahead[:end] if end else self.ahead


# Line endings when reading in text mode (i.e. universal newlines).
_TEXT_EOL = re.compile(rb"\r\n|\r|\n")


def _text_line_end(data: bytes, start: int, eof: bool) -> int | None:
    """Return where the line starting at ``start`` ends, as read in text mode.

    This is the offset just after its line ending. Returns None if that is not known
    from ``data``, the bytes at the start of a file (which reach the end of the file
    if ``eof``).
    """
    match = _TEXT_EOL.search(data, start)
    if match is None:
        return len(data) if eof else None
    # If this is a "\r" at the end of the data, it may be the start of a "\r\n". But
    # then nothing after it is known, so the head can't end here.
    return match.end()


def _head_end(data: bytes, eof: bool) -> int | None:
    """Return where the head of a file ends, given the bytes at its start.

    The head is everything up to the first comment line, that line, and the front
    matter of the data line after it (up to and including " spectrum ", or the whole
    line if there is no such marker). Returns None if more of the file is needed to
    find its end.
    """
    start = 0
    while (end := _text_line_end(data, start, eof)) is not None:
        if end == start:
            # The end of the file, without a comment line.
            return end
        if data.startswith(b"#", start):
            break
        start = end
    else:
        return None

    data_end = _text_line_end(data, end, eof)
    sep = _SPECTRUM_SEP.encode("ascii")
    marker = data.find(sep, end, len(data) if data_end is None else data_end)
    if marker >= 0:
        return marker + len(sep)
    return data_end


def _read_head(reader: _Reader) -> list[str]:
    """Read the head of an ACQ file, as lines of text.

    The head is the header, and the first comment line and front matter of the data
    line after it (see :func:`_head_end`). It is read in binary, but the lines are as
    they would be when reading the file in text mode.

    ``reader`` must be at the start of the file, and is left there, keeping the bytes
    it has read.
    """
    size = _HEAD_READ_SIZE
    while (end := _head_end(reader.peek(size), reader.eof)) is None:
        size *= 2
    return io.TextIOWrapper(io.BytesIO(reader.ahead[:end])).readlines()


def _iter_cycles(
    fl,
    fastspec: bool,
    read_spectrum: bool = True,
    progress: bool = True,
    leave_progress: bool = True,
    desc: str = "",
):
    """Iterate over the complete switch cycles in an open (binary-mode) ACQ file.

    Incomplete cycles and malformed lines are skipped (with a warning for the latter).

    Yields
    ------
    datas : tuple of DataEntry
        The three entries (swpos 0, 1, 2) of a complete cycle. Their spectra are
        None unless ``read_spectrum`` is True.
    offsets : tuple of int
        The byte offsets of the start of each of the three data lines in the file.
    """
    # First find the first swpos=0 line.
    for line in fl:
        if line.startswith(b"#"):
            cline = CommentLine.read(line.decode("ascii"), fastspec=fastspec)
            if cline.swpos == 0:
                break
    else:
        # No swpos=0 comment line, so there are no complete cycles.
        return

    offset = fl.tell()
    try:
        dline = next(fl)
    except StopIteration:
        # The file ends straight after the first swpos=0 comment line.
        return
    data = DataLine.read(dline.decode("ascii"), read_spectrum=read_spectrum)
    data = DataEntry(comment=cline, data=data)

    datas = (data,)
    offsets = (offset,)
    for line in tqdm.tqdm(
        fl,
        disable=not progress,
        desc=desc,
        unit="lines",
        leave=leave_progress,
    ):
        if line.startswith(b"\x00"):
            # The rest of the file is NUL-padded (e.g. an interrupted write).
            _warn_nul_padded(fl.name)
            break

        try:
            cline = CommentLine.read(line.decode("ascii"), fastspec=fastspec)
        except ACQError as e:
            warnings.warn(str(e), stacklevel=1)
            datas = ()
            continue

        offset = fl.tell()
        try:
            dline = next(fl).decode("ascii")
            data = DataLine.read(dline, read_spectrum=read_spectrum)
        except StopIteration:
            # We reached the end of the file.
            break
        except ACQLineError as e:
            # Something was bad in this line. Remove this iteration from the data
            # But try to keep going in the file.
            warnings.warn(str(e), stacklevel=1)
            datas = ()
            continue

        try:
            data = DataEntry(comment=cline, data=data)
            # Without the spectrum, DataEntry can't validate its length, so do it here.
            if not read_spectrum and _encoded_length(dline) != cline.nspec:
                raise ACQLineError("nspec and length of spectrum do not match")
        except ACQLineError as e:
            warnings.warn(str(e), stacklevel=1)
            datas = ()
            continue

        if data.comment.swpos == len(datas):
            # Add this to the cycle
            datas += (data,)
            offsets = (*offsets[: len(datas) - 1], offset)
        elif data.comment.swpos == 0:
            # Discard the previous cycle and start again -- it was incomplete.
            datas = (data,)
            offsets = (offset,)
        else:
            # Discard everything and keep going.
            datas = ()
            continue

        if data.comment.swpos == 2:
            # We have a full cycle
            yield datas, offsets
            datas = ()


def decode_file(
    fname: str | Path,
    progress: bool = True,
    leave_progress: bool = True,
):
    """
    Parse and decode an ACQ file, optionally writing it to a new format.

    Parameters
    ----------
    fname : str or Path
        filename of the ACQ file to read.

    progress: bool, optional
        Whether to display a progress bar for the read.
    meta: bool, optional
        Whether to output metadata for the read. Deprecated, will be set to True in a
        future version.
    leave_progress : bool, optional
        Whether to leave the progress bar (if one is used) on the screen when done.
        Useful to set to False if reading multiple files.
    """
    fname = Path(fname)
    p0, p1, p2 = [], [], []

    with fname.open("rb") as fl:
        anc = Ancillary._from_reader(fname, _Reader(fl))
        fastspec = "data_drops" in anc.meta
        fl.seek(0)
        for datas, _ in _iter_cycles(
            fl,
            fastspec=fastspec,
            progress=progress,
            leave_progress=leave_progress,
            desc=f"Reading {fname.name}",
        ):
            p0.append(datas[0].data.spectrum)
            p1.append(datas[1].data.spectrum)
            p2.append(datas[2].data.spectrum)
            anc.append(datas)

    anc.complete()
    p0 = np.array(p0)
    p1 = np.array(p1)
    p2 = np.array(p2)

    # Get Q ratio -- need to get this to have compatibility of output with new
    # fastspec default output (HDF5).
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=RuntimeWarning)
        q = (p0 - p1) / (p2 - p1)

    return q.T, [p0.T, p1.T, p2.T], anc


def _index_file(
    fname: str | Path,
    progress: bool = False,
    leave_progress: bool = True,
) -> tuple[Ancillary, np.ndarray]:
    """Read the ancillary data of an ACQ file without decoding any spectra.

    Returns
    -------
    anc : Ancillary
        The (completed) ancillary data, with one row per complete switch cycle, exactly
        as returned by :func:`decode_file`.
    offsets : np.ndarray
        An integer array of shape ``(ncycles, 3)`` giving the byte offset of the data
        line for each switch position of each cycle. Pass (a subset of) these to
        :func:`_read_spectra` to decode the spectra.
    """
    fname = Path(fname)
    offsets = []

    with fname.open("rb") as fl:
        anc = Ancillary._from_reader(fname, _Reader(fl))
        fastspec = "data_drops" in anc.meta
        fl.seek(0)
        for datas, offs in _iter_cycles(
            fl,
            fastspec=fastspec,
            read_spectrum=False,
            progress=progress,
            leave_progress=leave_progress,
            desc=f"Indexing {fname.name}",
        ):
            anc.append(datas)
            offsets.append(offs)

    anc.complete()
    return anc, np.array(offsets, dtype=np.int64).reshape((-1, 3))


def _read_spectra(
    fname: str | Path,
    offsets: np.ndarray,
    nchannels: int,
    channels: slice = slice(None),
) -> np.ndarray:
    """Decode the spectra at given byte offsets of an ACQ file.

    Parameters
    ----------
    fname : str or Path
        The ACQ file.
    offsets : np.ndarray
        Integer array of shape ``(ntimes, 3)`` of data-line offsets, as returned by
        :func:`_index_file`.
    nchannels : int
        The total number of channels in each spectrum.
    channels : slice
        A contiguous (unit-step) range of channels to decode.

    Returns
    -------
    np.ndarray
        The decoded spectra, with shape ``(3, ntimes, nchannels_selected)``.
    """
    start, stop, step = channels.indices(nchannels)
    if step != 1:
        raise ValueError("channels must be a contiguous slice")

    out = np.zeros((offsets.shape[1], offsets.shape[0], max(stop - start, 0)))
    sep = _SPECTRUM_SEP.encode("ascii")

    with Path(fname).open("rb") as fl:
        for i, offs in enumerate(offsets):
            for j, off in enumerate(offs):
                fl.seek(off)
                line = fl.readline()
                enc = line[line.index(sep) + len(sep) :].lstrip()
                _decode_into(enc[4 * start : 4 * stop], out[j, i])

    return out


def _complete_cycles(
    entries: Iterable[DataEntry | None],
) -> Iterator[tuple[DataEntry, DataEntry, DataEntry]]:
    """Group a stream of entries into complete (swpos 0, 1, 2) cycles.

    A ``None`` in the stream marks a bad entry, and discards the cycle in progress.
    This is the same grouping as in :func:`_iter_cycles`.
    """
    datas = ()
    for data in entries:
        if data is None:
            datas = ()
            continue

        if data.comment.swpos == len(datas):
            # Add this to the cycle
            datas += (data,)
        elif data.comment.swpos == 0:
            # Discard the previous cycle and start again -- it was incomplete.
            datas = (data,)
        else:
            # Discard everything and keep going.
            datas = ()
            continue

        if data.comment.swpos == 2:
            # We have a full cycle
            yield datas
            datas = ()


# Bytes read at the start of each entry: enough for the comment line and the front
# matter of the data line (together ~160 bytes). If that's not enough, up to
# _MAX_HEAD_SIZE bytes are read.
_HEAD_SIZE = 256
_MAX_HEAD_SIZE = 4096


def _make_entry(cline: CommentLine, dline: DataLine, nchannels: int) -> DataEntry:
    """Create an entry whose spectrum was not decoded, checking it as DataEntry does.

    ``nchannels`` is the number of channels that decoding the spectrum would give.
    """
    entry = DataEntry(comment=cline, data=dline)
    if nchannels != cline.nspec:
        raise ACQLineError("nspec and length of spectrum do not match")
    return entry


def _read_entry(cline: CommentLine, line: str) -> DataEntry:
    """Create an entry from a full data line, without decoding the spectrum."""
    dline = DataLine.read(line, read_spectrum=False)
    return _make_entry(cline, dline, _encoded_length(line))


def _split_head(chunk: bytes) -> tuple[bytes, bytes, int] | None:
    """Split the start of an entry into its comment line and data-line front matter.

    Also returns the offset of the start of the encoded spectrum. Returns None if the
    chunk does not contain all of these.
    """
    comment_end = chunk.find(b"\n") + 1
    if not comment_end:
        return None

    sep = _SPECTRUM_SEP.encode("ascii")
    marker = chunk.find(sep, comment_end)
    if marker < 0 or b"\n" in chunk[comment_end:marker]:
        return None

    spec_start = len(chunk) - len(chunk[marker + len(sep) :].lstrip(b" "))
    if spec_start == len(chunk) or chunk[spec_start : spec_start + 1].isspace():
        return None

    return chunk[:comment_end], chunk[comment_end:marker], spec_start


def _line_ending(tail: bytes) -> int:
    """Length of the line ending at the start of ``tail``, or 0 if there is none.

    The line ending only counts if it is followed by a comment line or the end of the
    file, as it must be at the end of a data line.
    """
    for eol in (b"\n", b"\r\n"):
        if tail.startswith(eol) and tail[len(eol) : len(eol) + 1] in (b"#", b""):
            return len(eol)
    return 0


def _entry_head(reader: _Reader) -> tuple[bytes, tuple[bytes, int] | None]:
    """Read the head of the entry at the reader's position.

    Returns its comment line and, if they were found in a small read, the front matter
    of its data line and the offset of the start of its spectrum.
    """
    head = _split_head(reader.peek(_HEAD_SIZE))
    if head is None:
        head = _split_head(reader.peek(_MAX_HEAD_SIZE))
    if head is None:
        return reader.readline(), None
    comment, front, spec_start = head
    return comment, (front, spec_start)


def _skip_spectrum(reader: _Reader, spec_start: int, nspec: int) -> bool:
    """Move past the entry at the reader's position, if its data line ends as expected.

    Since the encoded spectrum of a good data line has exactly 4*nspec characters, we
    seek straight to where the line should end, and check that it does (reading the
    head of the next entry at the same time). A data line that is too short would be
    missed only if the lines after it happen to end exactly where it should have
    ended, which needs corruption spanning exactly whole lines.
    """
    end = reader.pos + spec_start + 4 * nspec
    tail, eof = reader.read(end, _HEAD_SIZE + 2)
    if eol := _line_ending(tail):
        reader.move(end + eol, tail[eol:], eof)
        return True
    return False


def _next_entry(
    reader: _Reader,
    cline: CommentLine,
    comment: bytes,
    head: tuple[bytes, int] | None,
) -> DataEntry | None:
    """Read the entry at the reader's position, without its spectrum, and move past it.

    ``comment`` and ``head`` are from :func:`_entry_head`, and ``cline`` is the parsed
    comment line. Returns None at the end of the file, and raises ACQLineError if the
    entry is bad.
    """
    if head is not None:
        front, spec_start = head
        dline = DataLine.read(f"{front.decode('ascii')} spectrum ", read_spectrum=False)
        if _skip_spectrum(reader, spec_start, cline.nspec):
            return DataEntry(comment=cline, data=dline)

    # The data line is not the expected length: read all of it.
    reader.skip(len(comment))
    line = reader.readline()
    reader.skip(len(line))
    if not line:
        # We reached the end of the file.
        return None
    return _read_entry(cline, line.decode("ascii"))


def _iter_entries_without_spectra(
    fname: Path, reader: _Reader, fastspec: bool
) -> Iterator[DataEntry | None]:
    """Yield the entries of an ACQ file without reading their spectra.

    This yields the entries (and warnings) that :func:`_iter_cycles` reads, but with
    no spectra, reading only the head of each entry wherever possible (see
    :func:`_skip_spectrum`). If an entry does not have the expected form, we fall back
    to reading it line by line. ``None`` is yielded for a bad entry.

    ``reader`` must be at the start of the file, and ``fname`` is its name (for
    warnings).
    """
    # As in _iter_cycles, start at the first swpos=0 entry. Before that, every comment
    # line is parsed (so a bad one is an error), and other lines are skipped.
    while True:
        if not reader.peek(_HEAD_SIZE).startswith(b"#"):
            line = reader.readline()
            if not line:
                return
            reader.skip(len(line))
            continue

        comment, head = _entry_head(reader)
        cline = CommentLine.read(comment.decode("ascii"), fastspec=fastspec)
        if cline.swpos == 0:
            break

        # Skip the data line, unless it could be read as a comment line.
        if (
            head is None
            or head[0].startswith(b"#")
            or not _skip_spectrum(reader, head[1], cline.nspec)
        ):
            reader.skip(len(comment))

    # As in _iter_cycles, a bad first entry is an error.
    entry = _next_entry(reader, cline, comment, head)
    if entry is None:
        return
    yield entry

    while True:
        if reader.peek(_HEAD_SIZE).startswith(b"\x00"):
            # The rest of the file is NUL-padded (e.g. an interrupted write).
            _warn_nul_padded(fname)
            return

        comment, head = _entry_head(reader)
        if not comment:
            return

        try:
            cline = CommentLine.read(comment.decode("ascii"), fastspec=fastspec)
        except ACQError as e:
            # As in decode_file, skip just this line, and try the next as a comment.
            warnings.warn(str(e), stacklevel=1)
            yield None
            reader.skip(len(comment))
            continue

        try:
            entry = _next_entry(reader, cline, comment, head)
        except ACQLineError as e:
            warnings.warn(str(e), stacklevel=1)
            entry = None
        else:
            if entry is None:
                # We reached the end of the file.
                return
        yield entry


def read_metadata(fname: str | Path) -> tuple[dict, dict[str, np.ndarray]]:
    """Read the metadata and per-cycle ancillary data of an ACQ file.

    This gives the same metadata and ancillary data as :func:`decode_file` (i.e.
    ``anc.meta`` and ``anc.data`` of the :class:`Ancillary` it returns), including
    dropping incomplete or bad cycles in the same way. However, it neither decodes
    nor reads the spectra: it seeks over them, so it is much faster, and reads only a
    small fraction of the file.

    Parameters
    ----------
    fname : str or Path
        filename of the ACQ file to read.

    Returns
    -------
    meta : dict
        The metadata of the file, from its header and first entry.
    ancillary : dict
        Arrays of the ancillary data of each complete cycle, with shape
        ``(ncycles, 3)``: "times", "adcmax", "adcmin" and (for fastspec files)
        "data_drops".
    """
    fname = Path(fname)

    # Unbuffered, so that each read reads only what is asked for.
    with fname.open("rb", buffering=0) as fl:
        reader = _Reader(fl)
        anc = Ancillary._from_reader(fname, reader)
        fastspec = "data_drops" in anc.meta
        entries = _iter_entries_without_spectra(fname, reader, fastspec)
        for datas in _complete_cycles(entries):
            anc.append(datas)

    anc.complete()
    return anc.meta, anc.data


def encode(
    filename: str | Path,
    p: list[np.ndarray],
    meta: dict,
    ancillary: dict[str, np.ndarray],
):
    """
    Encode raw powers and ancillary data as an ACQ file.

    Parameters
    ----------
    filename : path
        Path to output file to write.
    p : list of ndarray
        List of three ndarrays, one for each switch. Each array should be 2D, shape
        Ntimes x Nfreq.
    meta : dict
        Dictionary of metadata associated with the file, to be written to the header
        and preambles.
    ancillary : structured array
        Time-dependent ancillary information, such as times and adcmin/adcmax.
    """
    p = np.array(p)

    # data_drops is optional in the .h5 file because we don't really need it.
    data_drops = ancillary.get("data_drops", 0)

    time_has_3dim = len(ancillary["times"].shape) == 2

    temperature = int(meta["temperature"])
    nblk = meta["nblk"]
    nfreq = meta["nfreq"]
    meta = {k: v for k, v in meta.items() if k not in ["temperature", "nblk", "nfreq"]}

    with Path(filename).open("w") as fl:
        # Write the header
        fl.writelines(f";--{k}: {v}\n" for k, v in meta.items())

        # Go through each time
        for i in range(len(p[0])):
            for switch, pp in enumerate(p[:, i]):
                dd = (
                    data_drops[i, switch]
                    if hasattr(data_drops, "__len__")
                    else data_drops
                )
                fl.write(
                    f"# swpos {switch} "
                    f"data_drops {int(dd):>4} "
                    f"adcmax  {ancillary['adcmax'][i, switch]:.5f} "
                    f"adcmin {ancillary['adcmin'][i, switch]:.5f} "
                    f"temp  {temperature} C "
                    f"nblk {nblk} "
                    f"nspec {nfreq}\n"
                )

                time = ancillary["times"][i]

                if time_has_3dim:
                    time = time[switch]

                if isinstance(time, np.bytes_):
                    time = time.decode()

                fl.write(
                    f"{time} {switch} {meta['freq_min']}  "
                    f"{meta['freq_res']}  {meta['freq_max']}  "
                    f"0.3 spectrum "
                )

                fl.write(_encode(pp))
                fl.write("\n")
