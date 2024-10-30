"""Perform tests of the _coordinates_alan module."""

from datetime import datetime

import numpy as np
from astropy import units as un
from astropy.coordinates import Galactic, SkyCoord

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


def test_galactic_to_radec():
    ra = 0.7
    dec = 0.0

    coord = SkyCoord(ra=ra * un.rad, dec=dec * un.rad)
    galframe = Galactic()
    gal = coord.transform_to(galframe)

    _ra, _dec = crda.galactic_to_radec(gal.b.deg, gal.l.deg)
    assert np.isclose(_dec, dec, atol=1e-2)
    assert np.isclose(_ra, ra, atol=1e-2)


def test_non_leap_year_first_day():
    secs = (365 + 365) * 86400 + 1
    year, day, hour, minute, sec = crda.toyrday(secs)
    assert year == 1972
    assert day == 1
    assert hour == 0
    assert minute == 0
    assert sec == 1


def test_datetime_to_secs():
    dt = datetime(year=1970, month=1, day=1, hour=0, minute=0, second=0)
    assert crda.datetime_tosecs(dt) == 0

    dt = datetime(year=1970, month=1, day=1, hour=1, minute=0, second=0)
    assert crda.datetime_tosecs(dt) == 3600

    dt = datetime(year=1970, month=1, day=2, hour=0, minute=0, second=0)
    assert crda.datetime_tosecs(dt) == 86400
