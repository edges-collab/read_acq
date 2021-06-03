import re
import warnings
from pathlib import Path
from typing import List
from .codec import _encode_line, _decode_line

import numpy as np
import tqdm
from . import writers

from edges_io.h5 import HDF5RawSpectrum


try:
    import h5py

    HAVE_H5PY = True
except ImportError:
    HAVE_H5PY = False


class Ancillary:
    """Represents ancillary data of a file"""

    header_char = ";"
    DTYPE = np.dtype(
        [
            ("adcmax", np.float32),
            ("adcmin", np.float32),
            ("times", "S17"),
        ]
    )

    _splits = re.compile(r"[\d\.:]+")

    def __init__(self, fname):
        self.fastspec_version = self._get_fastspec_version(fname)

        self.check_file(fname)

        self.size, nlines = self.get_ntimes(fname)

        self.meta = self.read_metadata(fname)

        self.meta["n_file_lines"] = nlines

        self.data = {
            "adcmax": np.zeros((self.size, 3), dtype=np.float32),
            "adcmin": np.zeros((self.size, 3), dtype=np.float32),
            "times": np.zeros(self.size, dtype="S17"),
        }
        if "data_drops" in self.meta:
            self.data["data_drops"] = np.zeros((self.size, 3), dtype=int)

        self._current_size = 0

    def check_file(self, fname):
        with open(fname) as fl:
            cline = ""
            while not cline.startswith("#"):
                cline = fl.readline()

        cline = self._splits.findall(cline)
        if int(cline[0]) != 0:
            raise OSError(
                "The format of the ACQ file is incorrect: should start with swpos = 0"
            )

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

            for line in fl.readlines():

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
        out = self._read_header(fname)

        with open(fname) as fl:
            cline = ""
            while not cline.startswith("#"):
                cline = fl.readline()
            dateline = fl.readline()

        clines = self._splits.findall(cline)
        dateline = self._splits.findall(dateline)

        out.update(
            {
                "data_drops"
                if "data_drops" in clines
                else "resolution": float(clines[1]),
                "temperature": float(clines[4]),
                "nblk": int(clines[5]),
                "nfreq": int(clines[6]),
                "freq_min": float(dateline[2]),
                "freq_max": float(dateline[4]),
                "freq_res": float(dateline[3]),
            }
        )
        return out

    @property
    def frequencies(self):
        return np.linspace(
            self.meta["freq_min"], self.meta["freq_max"], self.meta["nfreq"]
        )

    def get_ntimes(self, fname):
        # count lines
        ntimes = 0
        with open(fname) as fl:
            nlines = 0
            for line in fl.readlines():
                if line[0] not in "#*;" and line:
                    ntimes += 1
                nlines += 1
        return ntimes // 3, nlines

    @property
    def times(self):
        """The times associated with the ancillary data as datetime objects"""
        raise NotImplementedError("We haven't done this yet.")

    def _add_to_self(self, add_to_self=None):
        if add_to_self is None:
            add_to_self = self._current_size < self.size
        elif add_to_self and self._current_size == self.size:
            raise ValueError("You cannot add any more lines to this object")
        return add_to_self

    def parse_cline(self, line, add_to_self=None):
        add_to_self = self._add_to_self(add_to_self)
        line = self._splits.findall(line)

        if add_to_self:
            swpos = int(line[0])
            if "data_drops" in self.data:
                self.data["data_drops"][self._current_size, swpos] = line[1]
            self.data["adcmax"][self._current_size, swpos] = line[2]
            self.data["adcmin"][self._current_size, swpos] = line[3]

        else:
            out = {
                "adcmax": float(line[2]),
                "adcmin": float(line[3]),
            }
            if "data_drops" in self.data:
                out["data_drops"] = int(line[1])
            return out

    def parse_specline(self, line, add_to_self=None):
        """Parse an ancillary data line of a file, and return a dict"""
        add_to_self = self._add_to_self(add_to_self)

        line = self._splits.findall(line)

        if add_to_self:
            self.data["times"][self._current_size] = line[0]
            self._current_size += 1
        else:
            return [line[0]] + line[2:]


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
    n_times = anc.size

    p = [
        np.empty((n_times, anc.meta["nfreq"])),
        np.empty((n_times, anc.meta["nfreq"])),
        np.empty((n_times, anc.meta["nfreq"])),
    ]

    with open(fname) as fl:
        i = 0
        for line in tqdm.tqdm(
            fl.readlines(),
            disable=not progress,
            total=anc.meta["n_file_lines"],
            desc=f"Reading {Path(fname).name}",
            unit="lines",
            leave=leave_progress,
        ):
            switch_state = i % 3
            i_time = i // 3

            # Break when we hit the last "2" switch state (there may be a dangling 0 and 1)
            if i_time == n_times:
                break

            if not line.strip() or line[0] in "*;":  # Empty line.
                continue

            if line.startswith("#"):  # comment line
                anc.parse_cline(line)
                continue

            try:
                line_ancillary, data_string = line.split("spectrum")
            except ValueError:
                # Uh-oh, we probably have an incomplete file.
                raise Exception(
                    "Formatting of ACQ file is incorrect on line starting\n{}".format(
                        line[:50]
                    )
                )

            # Read the spectrum
            spec = _decode_line(data_string.lstrip())

            # Parse ancillary data
            if not switch_state:
                anc.parse_specline(line_ancillary)

            # Save current spectrum into p at the right switch state
            p[switch_state][i_time] = spec

            i += 1

    # Get Q ratio -- need to get this to have compatibility of output with new
    # fastspec default output (HDF5).
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=RuntimeWarning)
        Q = (p[0] - p[1]) / (p[2] - p[1])

    if meta:
        return Q.T, [pp.T for pp in p], anc
    else:
        return Q, [pp.T for pp in p]


def decode_files(files, *args, **kwargs):
    """Call :func:`decode_file` on a list of files."""
    for fl in tqdm.tqdm(
        files, disable=len(files) < 5, desc="Processing files", unit="files"
    ):
        decode_file(fl, progress=len(files) < 5, *args, **kwargs)


def convert_file(
    fname: [str, Path],
    outfile: [str, None, Path] = None,
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
        size=anc.size,
        **{f"p{i}": p[i] for i in range(3)},
    )
    return Q, p, anc


def encode(
    filename: [str, Path], p: List[np.ndarray], meta: dict, ancillary: np.ndarray
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

    with open(filename, "w") as fl:
        # Write the header
        for k, v in meta.items():
            fl.write(f";--{k}: {v}\n")

        # Go through each time
        for i in range(len(p[0])):
            for switch, pp in enumerate(p[:, i]):
                fl.write(
                    f"# swpos {switch} data_drops\t{ancillary['data_drops'][i, switch]} "
                    f"adcmax  {ancillary['adcmax'][i, switch]} "
                    f"adcmin {ancillary['adcmin'][i, switch]} "
                    f"temp  {meta['temp']} C "
                    f"nblk {meta['nblk']} "
                    f"nspec {meta['nfreq']}\n"
                )

                fl.write(
                    f"{ancillary['time'][i]} {switch}\t{meta['freq_min']}  "
                    f"{meta['freq_res']}  {meta['freq_max']}  "
                    f"0.3 spectrum "
                )

                fl.write(_encode_line(pp, meta["nblk"]))


def convert_h5(infile, outfile):
    """Convert a HDF5 file into an ACQ file.

    HDF5 file must be in the new format used by fastspec.
    """
    obj = HDF5RawSpectrum(infile)

    p = [obj["spectra"]["p0"], obj["spectra"]["p1"], obj["spectra"]["p2"]]
    meta = obj["meta"]
    ancillary = np.zeros(
        (len(p[0]), 3),
        dtype=[("adcmax", float), ("adcmin", float), ("data_drops", int)],
    )
    ancillary["adcmax"] = obj["time_ancillary"]["adcmax"]
    ancillary["adcmin"] = obj["time_ancillary"]["adcmin"]
    ancillary["data_drops"] = obj["time_ancillary"]["data_drops"]

    encode(outfile, p=p, meta=meta, ancillary=ancillary)
