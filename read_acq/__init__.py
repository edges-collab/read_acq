import ctypes
import glob
import os
import re
import warnings
from os import path
from pathlib import Path

import numpy as np
import tqdm
from . import writers

try:
    import h5py

    HAVE_H5PY = True
except ImportError:
    HAVE_H5PY = False

cdll = glob.glob(os.path.join(os.path.dirname(os.path.abspath(__file__)), "decode.*.so"))[0]
cdll = ctypes.CDLL(cdll)


def _decode_line(line):
    """
    Decode a pre-parsed line of an ACQ file.

    This is a very thin wrapper around the C `decode` function which does the same.

    :param line: the string line of the file, pre-parsed.
    :return: 1D-array of decoded values.
    """
    fnc = cdll.decode
    fnc.restype = ctypes.c_int
    fnc.argtypes = [
        ctypes.c_char_p,
        np.ctypeslib.ndpointer(np.float64),
    ]

    out = np.zeros(len(line) // 4)

    res = fnc(ctypes.c_char_p(line.encode("ascii")), np.ascontiguousarray(out))

    if res:
        raise Exception("C decoder exited with an error!")

    return out


PART_COEFFS = {
    'TR136_170': (1.03514e-3, 2.33825e-4, 7.92467e-8),
}


def convert_thermistor_ohm_to_kelvin(r, partname='TR136_170'):
    """
    Convert thermistor resistance to temperature

    Parameters
    ----------
    r : float
        Resistance
    partname: str, optional
        The name of the thermistor product. Must exist in `PART_COEFFS`

    Returns
    -------
    float: temperature in K.
    """
    a = PART_COEFFS[partname]
    return 1 / (a[0] + a[1] * np.log(r) + a[2] * np.log(r) ** 3)


def _get_tcal(fname, tcal=None):
    if tcal is None:
        dr = os.path.join(os.path.dirname(fname), path.pardir, "Resistance")
        try:
            tcal = glob.glob(os.path.join(dr, "*.csv"))[0]
        except IndexError:
            raise ValueError("tcal not given, and no relevant tcal file found at ../Restistance/")

    if type(tcal) == str:
        tcal = Path(tcal)

    if isinstance(tcal, Path):
        data = np.genfromtxt(tcal, usecols=3, delimiter=',')
        tcal = np.mean(data[len(data) // 20:])
        tcal = convert_thermistor_ohm_to_kelvin(tcal)

    return tcal


class Ancillary:
    """Represents ancillary data of a file"""
    header_char = ";"
    DTYPE = np.dtype([
        ('adcmax', np.float32),
        ('adcmin', np.float32),
        ("time", "S17"),
    ])

    _splits = re.compile(r"[\d\.:]+")

    def __init__(self, fname):
        self.check_file(fname)

        self.size, nlines = self.get_ntimes(fname)
        self.meta = self.read_metadata(fname)

        self.meta['n_file_lines'] = nlines

        self.data = np.zeros(self.size, dtype=self.DTYPE)
        self._current_size = 0

    def check_file(self, fname):
        with open(fname) as fl:
            cline = ''
            while not cline.startswith('#'):
                cline = fl.readline()

        cline = self._splits.findall(cline)
        if int(cline[0]) != 0:
            raise IOError("The format of the ACQ file is incorrect: should start with swpos = 0")

    def _read_header(self, fname):
        out = {}

        name_pattern = re.compile(r"[a-zA-Z_]+")
        type_order = [int, float, str]

        with open(fname) as fl:
            for line in fl.readlines():
                if not line.startswith(self.header_char):
                    break
                name, val = line.split(":")
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
            cline = ''
            while not cline.startswith("#"):
                cline = fl.readline()
            dateline = fl.readline()

        cline = self._splits.findall(cline)
        dateline = self._splits.findall(dateline)

        out.update(
            {
                "resolution": float(cline[1]),
                "temperature": float(cline[4]),
                'nblk': int(cline[5]),
                'nfreq': int(cline[6]),
                'freq_min': float(dateline[2]),
                'freq_max': float(dateline[4]),
                'freq_res': float(dateline[3])
            }
        )
        return out

    @property
    def frequencies(self):
        return np.linspace(self.meta['freq_min'], self.meta['freq_max'], self.meta['nfreq'])

    def __getitem__(self, item):
        try:
            # treat it as an int
            return self.data[int(item)]
        except ValueError:
            # try it as a time
            return self.data[self.data['time'] == item]

    def get_ntimes(self, fname):
        # count lines
        ntimes = 0
        with open(fname, 'r') as fl:
            nlines = 0
            for line in fl.readlines():
                if line[0] not in '#*;' and line:
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
            self.data[self._current_size]['adcmax'] = line[2]
            self.data[self._current_size]['adcmin'] = line[3]
        else:
            return {"adcmax": float(line[2]), "admin": float(line[3])}

    def parse_specline(self, line, add_to_self=None):
        """Parse an ancillary data line of a file, and return a dict"""
        add_to_self = self._add_to_self(add_to_self)

        line = self._splits.findall(line)

        if add_to_self:
            self.data[self._current_size]['time'] = line[0]
            self._current_size += 1
        else:
            return [line[0]] + line[2:]


def decode_file(fname, tcal=400, tload=300, outfile=None, write_formats=None,
                progress=True, ):
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
    """
    tcal = _get_tcal(fname, tcal)

    for fmt in write_formats:
        if fmt not in writers._WRITERS:
            raise ValueError("The format {} does not have an associated writer.".format(fmt))

    anc = Ancillary(fname)
    ntimes = anc.size

    p = [np.empty((ntimes, anc.meta['nfreq']))] * 3

    i = 0
    with open(fname, 'r') as fl:
        for line in tqdm.tqdm(
                fl.readlines(), disable=not progress, total=anc.meta['n_file_lines'],
                desc="Reading {}".format(fname), unit='lines'
        ):
            switch_state = i % 3
            i_time = i // 3

            # Break when we hit the last "2" switch state (there may be a dangling 0 and 1)
            if i_time == ntimes:
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
                raise Exception('Formatting of ACQ file is incorrect on line starting\n{}'.format(line[:50]))

            # Read the spectrum
            spec = _decode_line(data_string.lstrip())

            # Parse ancillary data
            if not switch_state:
                anc.parse_specline(line_ancillary)

            # Check the switch state and act accordingly
            p[switch_state][i_time] = spec

            # Move on unless switch is 2.
            if switch_state < 2:
                continue

            # Check for consistency of data.
            if not len(set([len(xx) for xx in p])) == 1:
                # p are different lengths!
                raise Exception('Error: Failed to calibrate.')

            i_time += 1
            i += 1

    # Get antenna temperature
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=RuntimeWarning)
        ant_temp = (p[0] - p[1]) / (p[2] - p[1]) * tcal + tload

    if write_formats is None:
        write_formats = ['mat']

    for fmt in write_formats:
        getattr(writers, '_write_%s' % fmt)(
            outfile=outfile or fname,
            ancillary=anc.meta,
            ant_temp=ant_temp.T,
            time_data=anc.data,
            freqs=anc.frequencies,
            **{'p{}'.format(i): p[i].T for i in range(3)},
        )

    return ant_temp, p


def decode_files(files, *args, **kwargs):
    """Call :func:`decode_file` on a list of files."""
    for fl in tqdm.tqdm(files, disable=len(files) < 5, desc="Processing files", unit='files'):
        decode_file(fl, progress=len(files) < 5, *args, **kwargs)
