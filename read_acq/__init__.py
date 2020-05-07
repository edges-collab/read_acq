import re
import warnings
from pathlib import Path
from typing import List
from .codec import _encode_line, _decode_line

import numpy as np
import tqdm
from . import writers
from pkg_resources import get_distribution, DistributionNotFound

try:
    __version__ = get_distribution(__name__).version
except DistributionNotFound:
    # package is not installed
    pass

try:
    import h5py

    HAVE_H5PY = True
except ImportError:
    HAVE_H5PY = False


class Ancillary:
    """Represents ancillary data of a file"""

    header_char = ";"
    DTYPE = np.dtype([("adcmax", np.float32), ("adcmin", np.float32), ("time", "S17"),])

    _splits = re.compile(r"[\d\.:]+")

    def __init__(self, fname):
        self.fastspec_version = self._get_fastspec_version(fname)

        self.check_file(fname)

        self.size, nlines = self.get_ntimes(fname)

        self.meta = self.read_metadata(fname)

        self.meta["n_file_lines"] = nlines

        self.data = np.zeros(self.size, dtype=self.DTYPE)
        self._current_size = 0

    def check_file(self, fname):
        with open(fname) as fl:
            cline = ""
            while not cline.startswith("#"):
                cline = fl.readline()

        cline = self._splits.findall(cline)
        if int(cline[0]) != 0:
            raise IOError(
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
        type_order = [int, float, str]

        with open(fname) as fl:

            for line in fl.readlines():

                if not line.startswith(self.header_char):
                    break
                if line.startswith("; FASTSPEC"):
                    continue

                try:
                    name, val = line.split(":")
                except ValueError:
                    print(fname, line)
                    raise
                name = name_pattern.findall(name)[0]
                for tp in type_order:
                    try:
                        out[name] = tp(val.split()[0])
                        break
                    except (ValueError, IndexError):
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

    def __getitem__(self, item):
        try:
            # treat it as an int
            return self.data[int(item)]
        except ValueError:
            # try it as a time
            return self.data[self.data["time"] == item]

    def get_ntimes(self, fname):
        # count lines
        ntimes = 0
        with open(fname, "r") as fl:
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
            self.data[self._current_size]["adcmax"] = line[2]
            self.data[self._current_size]["adcmin"] = line[3]
        else:
            return {"adcmax": float(line[2]), "admin": float(line[3])}

    def parse_specline(self, line, add_to_self=None):
        """Parse an ancillary data line of a file, and return a dict"""
        add_to_self = self._add_to_self(add_to_self)

        line = self._splits.findall(line)

        if add_to_self:
            self.data[self._current_size]["time"] = line[0]
            self._current_size += 1
        else:
            return [line[0]] + line[2:]


def decode_file(
    fname,
    outfile=None,
    write_formats=None,
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
    tcal : float or str/Path, optional
        The calibration temperature. If given as a str or Path, will attempt to read
        the given CSV file and take an average cal temperature. By default, looks in
        `../Resistance/` for this file.
    tload: float, optional
        The load temperature
    nchannels: int, optional
        Number of channels in the data.
    outfile: str, optional
        Filename to which to write the data. Only used if `write_formats` is not empty.
    write_formats: list of str, optional
        Can contain any of 'mat' and 'npz'.
    progress: bool, optional
        Whether to display a progress bar for the read.
    meta: bool, optional
        Whether to output metadata for the read.
    """
    for fmt in write_formats:
        if fmt not in writers._WRITERS:
            raise ValueError(
                "The format {} does not have an associated writer.".format(fmt)
            )

    anc = Ancillary(fname)
    n_times = anc.size

    p = [
        np.empty((n_times, anc.meta["nfreq"])),
        np.empty((n_times, anc.meta["nfreq"])),
        np.empty((n_times, anc.meta["nfreq"])),
    ]

    i = 0
    with open(fname, "r") as fl:
        for line in tqdm.tqdm(
            fl.readlines(),
            disable=not progress,
            total=anc.meta["n_file_lines"],
            desc="Reading {}".format(Path(fname).name),
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

    # Get antenna temperature
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=RuntimeWarning)
        Q = (p[0] - p[1]) / (p[2] - p[1])

    if write_formats is None:
        write_formats = ["mat"]

    for fmt in write_formats:
        getattr(writers, "_write_%s" % fmt)(
            outfile=outfile or fname,
            ancillary=anc.meta,
            Qratio=Q.T,
            time_data=anc.data,
            freqs=anc.frequencies,
            fastspec_version=anc.fastspec_version,
            size=anc.size,
            **{"p{}".format(i): p[i].T for i in range(3)},
        )

    if meta:
        return Q, p, anc
    else:
        return Q, p


def decode_files(files, *args, **kwargs):
    """Call :func:`decode_file` on a list of files."""
    for fl in tqdm.tqdm(
        files, disable=len(files) < 5, desc="Processing files", unit="files"
    ):
        decode_file(fl, progress=len(files) < 5, *args, **kwargs)


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
                    f"# swpos {switch} data_drops\t{meta['data_drops']} "
                    f"adcmax  {ancillary['adcmax'][i]} "
                    f"adcmin {ancillary['adcmin'][i]} "
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
