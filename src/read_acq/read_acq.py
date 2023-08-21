"""Functions and classes for reading and writing the .acq format."""
from __future__ import annotations

import attrs
import numpy as np
import re
import tqdm
import warnings
from pathlib import Path

from . import writers
from .codec import _decode_line, _encode_line


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
            raise ACQLineError(f"Could not parse line: '{line}' -- probably incomplete")

        match = re.match(cls.regex, front)
        if match is None:
            raise ACQError(f"Could not parse line front-matter: {front}")

        if read_spectrum:
            spec = _decode_line(back.lstrip())
        else:
            spec = None

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
    DTYPES = {
        "adcmax": np.float32,
        "adcmin": np.float32,
        "times": "S17",
    }

    def __init__(self, fname):
        self.fastspec_version = self._get_fastspec_version(fname)
        self.meta = self.read_metadata(fname)

        self.data = {
            "adcmax": [],  # np.zeros((self.size, 3), dtype=np.float32),
            "adcmin": [],  # np.zeros((self.size, 3), dtype=np.float32),
            "times": [],  # np.zeros((self.size, 3), dtype="S17"),
        }
        if "data_drops" in self.meta:
            self.data["data_drops"] = []  # np.zeros((self.size, 3), dtype=int)

    def _get_fastspec_version(self, fname):
        with open(fname) as fl:
            first_line = fl.readline()
        if first_line.startswith(self.header_char):
            return first_line.split("FASTSPEC")[-1]
        else:
            return None

    def _read_header(self, fname):
        out = {}

        name_pattern = re.compile(r"[a-zA-Z_]+")
        with open(fname) as fl:
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
                        warnings.warn(f"In file {fname}, item {line} has no value")
                        name = line.split(":")[0]
                        val = ""

                    name = name_pattern.findall(name)[0]

                for tp in type_order:
                    try:
                        out[name] = tp(val.split()[0])
                        break
                    except IndexError:
                        try:
                            out[name] = tp(val)
                        except ValueError:
                            pass
                    except ValueError:
                        pass

        return out

    def read_metadata(self, fname):
        """Read the metadata of the ACQ file."""
        out = self._read_header(fname)

        with open(fname) as fl:
            for line in fl:
                if line.startswith("#"):
                    comment = CommentLine.read(line)
                    data = DataLine.read(next(fl), read_spectrum=False)
                    break

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


def decode_file(
    fname,
    progress=True,
    meta=False,
    leave_progress=True,
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
        Whether to output metadata for the read.
    leave_progress : bool, optional
        Whether to leave the progress bar (if one is used) on the screen when done.
        Useful to set to False if reading multiple files.
    """
    anc = Ancillary(fname)

    fastspec = "data_drops" in anc.meta

    p0, p1, p2 = [], [], []

    with open(fname) as fl:
        # First find the first swpos=0 line.

        for line in fl:
            if line.startswith("#"):
                cline = CommentLine.read(line, fastspec=fastspec)
                if cline.swpos == 0:
                    break

        data = DataLine.read(next(fl))
        data = DataEntry(comment=cline, data=data)

        datas = (data,)
        for line in tqdm.tqdm(
            fl,
            disable=not progress,
            desc=f"Reading {Path(fname).name}",
            unit="lines",
            leave=leave_progress,
        ):
            cline = CommentLine.read(line, fastspec=fastspec)
            try:
                data = DataLine.read(next(fl))
            except StopIteration:
                # We reached the end of the file.
                break
            except ACQLineError as e:
                # Something was bad in this line. Remove this iteration from the data
                # But try to keep going in the file.
                warnings.warn(str(e))
                datas = ()
                continue

            try:
                data = DataEntry(comment=cline, data=data)
            except ACQLineError as e:
                warnings.warn(str(e))
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
                p0.append(datas[0].data.spectrum)
                p1.append(datas[1].data.spectrum)
                p2.append(datas[2].data.spectrum)
                anc.append(datas)
                datas = ()

    anc.complete()
    p0 = np.array(p0)
    p1 = np.array(p1)
    p2 = np.array(p2)

    # Get Q ratio -- need to get this to have compatibility of output with new
    # fastspec default output (HDF5).
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=RuntimeWarning)
        Q = (p0 - p1) / (p2 - p1)

    if meta:
        return Q.T, [p0.T, p1.T, p2.T], anc
    else:
        return Q, [p0.T, p1.T, p2.T]


def convert_file(
    fname: str | Path,
    outfile: str | None | Path = None,
    write_format: str = "h5",
    **kwargs,
):
    """Convert an ACQ file to a different format.

    Parameters
    ----------
    fname
        The filename of the ACQ file to convert.
    outfile
        The path to the output file. If None, save to a file of the same name (with
        different extension) as ``fname``.
    write_format
        The format to write to. Options are 'h5', 'mat' and 'npz'.

    Other Parameters
    ----------------
    All other parameters passed through to :func:`decode_file`.
    """
    kwargs["meta"] = True
    if write_format not in writers._WRITERS:
        raise ValueError(
            f"The format '{write_format}' does not have an associated writer."
        )

    Q, p, anc = decode_file(fname, **kwargs)

    getattr(writers, f"_write_{write_format}")(
        outfile=outfile or fname,
        ancillary=anc.meta,
        Qratio=Q,
        time_data=anc.data,
        freqs=anc.frequencies,
        fastspec_version=anc.fastspec_version,
        size=len(anc.data["adcmax"]),
        **{f"p{i}": p[i] for i in range(3)},
    )
    return Q, p, anc


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

    with open(filename, "w") as fl:
        # Write the header
        for k, v in meta.items():
            fl.write(f";--{k}: {v}\n")

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
                    f"data_drops {dd:>4} "
                    f"adcmax  {ancillary['adcmax'][i, switch]} "
                    f"adcmin {ancillary['adcmin'][i, switch]} "
                    f"temp  {meta['temp']} C "
                    f"nblk {meta['nblk']} "
                    f"nspec {meta['nfreq']}\n"
                )

                time = ancillary["times"][i]
                if time_has_3dim:
                    time = time[switch]

                fl.write(
                    f"{time} {switch} {meta['freq_min']}  "
                    f"{meta['freq_res']}  {meta['freq_max']}  "
                    f"0.3 spectrum "
                )

                fl.write(_encode_line(pp, meta["nblk"]))
                fl.write("\n")


def convert_h5(infile, outfile):
    """Convert a HDF5 file into an ACQ file.

    HDF5 file must be in the new format used by fastspec.
    """
    try:
        from edges_io.h5 import HDF5RawSpectrum
    except ImportError:  # pragma: no cover
        raise ImportError(
            "To use convert_h5, you need to install edges_io or do "
            "`pip install read_acq[h5]`"
        )

    obj = HDF5RawSpectrum(infile)

    p = [obj["spectra"]["p0"], obj["spectra"]["p1"], obj["spectra"]["p2"]]
    meta = obj["meta"]

    ancillary = {}
    ancillary["adcmax"] = obj["time_ancillary"]["adcmax"]
    ancillary["adcmin"] = obj["time_ancillary"]["adcmin"]
    ancillary["times"] = obj["time_ancillary"]["times"]

    encode(outfile, p=[pp.T for pp in p], meta=meta, ancillary=ancillary)
