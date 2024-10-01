"""Perform tests of the _coordinates_alan module."""

import numpy as np
from read_acq import _coordinates as crda


def test_alan_coordinates_toyrday():
    seconds_since_ny1970 = 1474043217.33333333
    year, day, hour, minute, sec = crda.toyrday(seconds_since_ny1970)
    assert year == 2016
    assert day == 260
    assert hour == 16


def test_ra_dec_at_lst():
    az, el = crda.radec_azel_from_lst(lst=0, ra=0, dec=0, lat=0)
    assert az == 0.0
    assert el == np.pi / 2
