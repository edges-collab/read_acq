"""Functions and classes for reading and writing the .acq format."""

from __future__ import annotations

import contextlib
import itertools
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
        self.fastspec_version = self._get_fastspec_version(fname)
        self.meta = self.read_metadata(fname)

        self.data = {
            "adcmax": [],  # np.zeros((self.size, 3), dtype=np.float32),
            "adcmin": [],  # np.zeros((self.size, 3), dtype=np.float32),
            "times": [],  # np.zeros((self.size, 3), dtype="S17"),
        }
        if "data_drops" in self.meta:
            self.data["data_drops"] = []  # np.zeros((self.size, 3), dtype=int)

    def _get_fastspec_version(self, fname: Path):
        with fname.open("r") as fl:
            first_line = fl.readline()
        if first_line.startswith(self.header_char):
            return first_line.split("FASTSPEC")[-1]
        else:
            return None

    def _read_header(self, fname: Path):
        out = {}

        name_pattern = re.compile(r"[a-zA-Z_]+")
        with fname.open("r") as fl:
            type_order = [int, float, str]

            for line in fl:
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
        out = self._read_header(fname)

        with fname.open("r") as fl:
            for line in fl:
                if line.startswith("#"):
                    comment = CommentLine.read(line)
                    data = DataLine.read(next(fl), read_spectrum=False)
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
    anc = Ancillary(fname)

    fastspec = "data_drops" in anc.meta
    p0, p1, p2 = [], [], []

    with fname.open("rb") as fl:
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
    anc = Ancillary(fname)

    fastspec = "data_drops" in anc.meta
    offsets = []

    with fname.open("rb") as fl:
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
# matter of the data line (together ~160 bytes).
_HEAD_SIZE = 512


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


def _iter_entries_without_spectra(
    fname: Path, fl: BinaryIO, fastspec: bool
) -> Iterator[DataEntry | None]:
    """Yield the entries of an ACQ file without reading their spectra.

    This yields the entries (and warnings) that :func:`_iter_cycles` reads, but with
    no spectra. Since the encoded spectrum of a good data line has exactly 4*nspec
    characters, we seek straight to where the line should end, and check that it does.
    If it doesn't, we fall back to reading the whole line. (A data line that is too
    short would be missed only if the lines after it happen to end exactly where it
    should have ended, which needs corruption spanning exactly whole lines.)

    ``fl`` must be a binary file positioned at the start of a comment line, and
    ``fname`` is its name (for warnings).
    """
    pos = fl.tell()
    while True:
        fl.seek(pos)
        chunk = fl.read(_HEAD_SIZE)
        if chunk.startswith(b"\x00"):
            # The rest of the file is NUL-padded (e.g. an interrupted write).
            _warn_nul_padded(fname)
            return

        head = _split_head(chunk)
        if head is None:
            # Not a complete, well-formed entry head (e.g. at the end of the file):
            # read line by line.
            fl.seek(pos)
            comment = fl.readline()
            if not comment:
                return
        else:
            comment, front, spec_start = head

        try:
            cline = CommentLine.read(comment.decode("ascii"), fastspec=fastspec)
        except ACQError as e:
            # As in decode_file, skip just this line, and try the next as a comment.
            warnings.warn(str(e), stacklevel=1)
            yield None
            pos += len(comment)
            continue

        if head is not None:
            dline = DataLine.read(f"{front.decode()} spectrum ", read_spectrum=False)
            end = pos + spec_start + 4 * cline.nspec
            fl.seek(end)
            if eol := _line_ending(fl.read(3)):
                try:
                    entry = DataEntry(comment=cline, data=dline)
                except ACQLineError as e:
                    warnings.warn(str(e), stacklevel=1)
                    entry = None
                yield entry
                pos = end + eol
                continue

        # The data line is not the expected length: read all of it.
        fl.seek(pos + len(comment))
        line = fl.readline()
        if not line:
            # We reached the end of the file.
            return
        try:
            entry = _read_entry(cline, line.decode("ascii"))
        except ACQLineError as e:
            warnings.warn(str(e), stacklevel=1)
            entry = None
        yield entry
        pos = fl.tell()


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
    anc = Ancillary(fname)
    fastspec = "data_drops" in anc.meta

    # A small buffer, so that each read at an entry does not read far beyond it.
    with fname.open("rb", buffering=2 * _HEAD_SIZE) as fl:
        # As in decode_file, start at the first swpos=0 entry.
        for line in fl:
            if line.startswith(b"#"):
                cline = CommentLine.read(line.decode("ascii"), fastspec=fastspec)
                if cline.swpos == 0:
                    break
        else:
            cline = None

        line = fl.readline() if cline is not None else b""
        if line:
            # As in decode_file, a bad first entry is an error.
            first = _read_entry(cline, line.decode("ascii"))
            entries = _iter_entries_without_spectra(fname, fl, fastspec)
            for datas in _complete_cycles(itertools.chain([first], entries)):
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
