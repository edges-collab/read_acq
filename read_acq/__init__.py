import copy
import ctypes
import glob
import os
import re
from os import path
from pathlib import Path

import numpy as np
import tqdm
from scipy import io

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


def _get_fname(infile=None, outfile=None):
    if infile and not outfile:
        outfile = path.splitext(infile)[0]
    elif not outfile and not infile:
        raise Exception("You need to provide either outfile or infile!")

    return outfile


def _write_npz(infile=None, outfile=None, **data):
    outfile = _get_fname(infile, outfile)
    np.savez(outfile, **data)


def _write_mat(infile=None, outfile=None, **data):
    outfile = _get_fname(infile, outfile)
    io.savemat(outfile + '.mat', data)


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
        if np.average(data)>55 and np.average(data)<45:
            tcal=350
        else:
            tcal = np.mean(data[len(data) // 20:])
            tcal = convert_thermistor_ohm_to_kelvin(tcal)
        #default should be 300 if data is between 45 and 55

    return tcal


def decode_file(fname, tcal=None, tload=300, nchannels=16384 * 2, outfile=None, write_formats=None,
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
    print("Tcal = ", tcal)
    print("Tload = ", tload)

    # count lines
    ntimes = 0
    with open(fname, 'r') as fl:
        for line in fl.readlines():
            if line[0] not in '#*' and line:
                switch_state = int(line.split()[1])
                if switch_state == 2:
                    ntimes += 1

    p = [np.empty((ntimes, nchannels)),
         np.empty((ntimes, nchannels)),
         np.empty((ntimes, nchannels))]

    anc = []

    i_time = 0
    splits = re.compile(r"[\w']+")
    with open(fname, 'r') as fl:
        for line in tqdm.tqdm(fl.readlines(), disable=not progress, total=3 * ntimes):
            # Break when we hit the last "2" switch state (there may be a dangling 0 and 1)
            if i_time == ntimes:
                break

            if not line or line[0] in '#*':  # Empty line.
                continue

            try:
                line_ancillary, data_string = line.split("spectrum")
            except:
                # Uh-oh, we probably have an incomplete file.
                raise Exception('Failed to parse ancillary data.')

            # Split ancillary into a list of ancillary data
            line_ancillary = splits.findall(line_ancillary)

            # Read the spectrum
            spec = _decode_line(data_string.lstrip())
            total_power = np.sum(spec)

            # Check the switch state and act accordingly
            switch_state = int(line_ancillary[5])
            p[switch_state][i_time] = spec

            if switch_state == 0:
                ancillary = copy.copy(line_ancillary)

            ancillary.append(total_power)

            # Move on unless switch is 2.
            if switch_state < 2:
                continue

            # Check for consistency of data.
            if not len(set([len(xx) for xx in p])) == 1:
                # p are different lengths!
                raise Exception('Error: Failed to calibrate.')

            anc.append(ancillary)

            i_time += 1

    if switch_state < 1:
        p[1] = p[1][:-1]
    if switch_state < 2:
        p[2] = p[2][:-1]

    # Get antenna temperature
    ant_temp = (p[0] - p[1]) / (p[2] - p[1]) * tcal + tload

    if write_formats is None:
        write_formats = ['mat']

    for fmt in write_formats:
        try:
            globals()['_write_%s' % fmt](
                infile=fname, outfile=outfile,
                ant_temp=ant_temp.T,
                **{'p{}'.format(i): p[i].T for i in range(3)},
            )
        except KeyError as e:
            print(e)
            raise NotImplementedError("the format {} is not supported".format(fmt))

    return ant_temp, p


def decode_files(files, *args, **kwargs):
    """Call :func:`decode_file` on a list of files."""
    for fl in tqdm.tqdm(files, disable=len(files) < 5):
        decode_file(fl, progress=len(files) < 5, *args, **kwargs)
