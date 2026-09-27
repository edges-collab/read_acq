"""Tests of reading only part of an ACQ file (selecting on read)."""

from pathlib import Path

import numpy as np
import pytest
from astropy import units as un
from astropy.time import Time
from pygsdata import GSData
from pygsdata.select import select_freqs, select_loads, select_lsts, select_times

from read_acq import decode_file, encode
from read_acq import read_acq as _ra
from read_acq.gsdata import fast_lst_setter, read_acq_to_gsdata
from read_acq.read_acq import ACQError, _index_file, _read_spectra

NTIMES = 30
NFREQ = 1024


@pytest.fixture(scope="module")
def long_acq(tmp_path_factory) -> Path:
    """Write an ACQ file with many switch cycles, spanning ~20 hours."""
    rng = np.random.default_rng(1234)
    start = Time("2023:070:00:00:00", format="yday", scale="utc")
    times = start + np.arange(NTIMES)[:, None] * 40 * un.min + [0, 13, 26] * un.s

    fname = tmp_path_factory.mktemp("acq") / "long.acq"
    encode(
        fname,
        p=list(rng.uniform(1e-3, 1, size=(3, NTIMES, NFREQ))),
        meta={
            "temperature": 25,
            "nblk": 2974,
            "nfreq": NFREQ,
            "freq_min": 0.0,
            "freq_res": 200.0 / NFREQ,
            "freq_max": 200.0,
        },
        ancillary={
            "times": np.array(
                [[t.strftime("%Y:%j:%H:%M:%S") for t in row] for row in times]
            ),
            "adcmax": rng.uniform(0, 1, size=(NTIMES, 3)),
            "adcmin": rng.uniform(-1, 0, size=(NTIMES, 3)),
        },
    )
    return fname


@pytest.fixture(scope="module", params=["sample", "pxspec", "long"])
def any_acq(request, sample_acq, long_acq) -> Path:
    if request.param == "long":
        return long_acq
    return sample_acq.parent / f"{request.param}.acq"


def _select_after_read(data: GSData, selectors: dict) -> GSData:
    """Apply selectors in the same order GSData.from_file does post-read."""
    if "freq_selector" in selectors:
        data = select_freqs(data, **selectors["freq_selector"])
    if "time_selector" in selectors:
        data = select_times(data, **selectors["time_selector"])
    if "lst_selector" in selectors:
        data = select_lsts(data, **selectors["lst_selector"])
    if "load_selector" in selectors:
        data = select_loads(data, **selectors["load_selector"])
    return data


def _assert_same(a: GSData, b: GSData):
    np.testing.assert_array_equal(a.data, b.data)
    np.testing.assert_array_equal(a.times.jd, b.times.jd)
    np.testing.assert_array_equal(a.time_ranges.jd, b.time_ranges.jd)
    np.testing.assert_allclose(a.lsts.hour, b.lsts.hour)
    np.testing.assert_allclose(a.lst_ranges.hour, b.lst_ranges.hour)
    np.testing.assert_array_equal(a.freqs, b.freqs)
    assert a.loads == b.loads
    assert a.name == b.name
    assert a.auxiliary_measurements.colnames == b.auxiliary_measurements.colnames
    for k in a.auxiliary_measurements.colnames:
        np.testing.assert_array_equal(
            a.auxiliary_measurements[k], b.auxiliary_measurements[k]
        )


def test_index_matches_decode(any_acq: Path):
    _, p, anc = decode_file(any_acq, progress=False)
    anc_idx, offsets = _index_file(any_acq)

    assert offsets.shape == (len(anc.data["times"]), 3)
    assert anc_idx.meta == anc.meta
    for k in anc.data:
        np.testing.assert_array_equal(anc_idx.data[k], anc.data[k])

    spec = _read_spectra(any_acq, offsets, anc.meta["nfreq"])
    np.testing.assert_array_equal(spec, np.transpose(p, (0, 2, 1)))

    part = _read_spectra(any_acq, offsets, anc.meta["nfreq"], slice(100, 900))
    np.testing.assert_array_equal(part, spec[..., 100:900])


def test_read_spectra_bad_slice(long_acq: Path):
    _, offsets = _index_file(long_acq)
    with pytest.raises(ValueError, match="contiguous"):
        _read_spectra(long_acq, offsets, NFREQ, slice(0, 100, 2))


def test_index_skips_truncated_lines(long_acq: Path, tmp_path: Path):
    lines = long_acq.read_text().splitlines(keepends=True)
    data_lines = [i for i, line in enumerate(lines) if " spectrum " in line]

    # Truncate the spectrum of switch 1 of the 5th cycle.
    bad = data_lines[5 * 3 + 1]
    lines[bad] = lines[bad][:-9] + "\n"
    fname = tmp_path / "truncated.acq"
    fname.write_text("".join(lines))

    with pytest.warns(UserWarning, match="nspec and length of spectrum do not match"):
        _, _, anc = decode_file(fname, progress=False)
    with pytest.warns(UserWarning, match="nspec and length of spectrum do not match"):
        anc_idx, offsets = _index_file(fname)

    assert len(anc.data["times"]) == NTIMES - 1
    np.testing.assert_array_equal(anc_idx.data["times"], anc.data["times"])
    assert len(offsets) == NTIMES - 1


T0 = Time("2023:070:00:00:00", format="yday", scale="utc")


@pytest.mark.parametrize(
    "selectors",
    [
        {"time_selector": {"indx": [2, 5, 7]}},
        {"time_selector": {"time_range": (T0 + 3 * un.hour, T0 + 9 * un.hour)}},
        {"lst_selector": {"lst_range": (6, 12)}},
        {"lst_selector": {"lst_range": (22, 2)}},
        {"lst_selector": {"lst_range": (6, 12), "gha": True}},
        {"freq_selector": {"freq_range": (50 * un.MHz, 100 * un.MHz)}},
        {"freq_selector": {"indx": [10, 20, 500]}},
        {"load_selector": {"loads": ("ant",)}},
        {
            "freq_selector": {"freq_range": (60 * un.MHz, 160 * un.MHz)},
            "time_selector": {"time_range": (T0 + 2 * un.hour, T0 + 15 * un.hour)},
            "lst_selector": {"indx": [0, 3, 4]},
            "load_selector": {"loads": ("ant", "internal_load")},
        },
    ],
)
@pytest.mark.parametrize("lst_setter", [None, fast_lst_setter])
def test_select_on_read_matches_select_after(long_acq, selectors, lst_setter):
    full = read_acq_to_gsdata(long_acq, lst_setter=lst_setter)
    expected = _select_after_read(full, selectors)
    got = read_acq_to_gsdata(long_acq, lst_setter=lst_setter, selectors=selectors)

    assert 0 < got.ntimes <= full.ntimes
    _assert_same(got, expected)


def test_multifile_indx_spans_files(long_acq):
    selectors = {
        "time_selector": {"indx": [NTIMES - 2, NTIMES - 1, NTIMES, NTIMES + 3]}
    }
    full = read_acq_to_gsdata([long_acq, long_acq])
    got = read_acq_to_gsdata([long_acq, long_acq], selectors=selectors)
    _assert_same(got, _select_after_read(full, selectors))


def test_multifile_skips_empty_selection_in_one_file(long_acq, sample_acq):
    # The sample file is from a different day, so this time range excludes it.
    selectors = {"time_selector": {"time_range": (T0, T0 + 5 * un.hour)}}
    got = read_acq_to_gsdata([long_acq, sample_acq], selectors=selectors)
    ref = read_acq_to_gsdata(long_acq, selectors=selectors)
    np.testing.assert_array_equal(got.data, ref.data)


def test_from_file_decodes_only_selection(long_acq, monkeypatch):
    calls = []
    decode_into = _ra._decode_into

    def spy(line, out):
        calls.append(len(line))
        decode_into(line, out)

    monkeypatch.setattr(_ra, "_decode_into", spy)

    gsd = GSData.from_file(
        long_acq,
        selectors={
            "time_selector": {"indx": [1, 2, 3]},
            "freq_selector": {"freq_range": (50 * un.MHz, 100 * un.MHz)},
        },
    )
    assert gsd.ntimes == 3
    assert len(calls) == 3 * 3
    assert all(n == 4 * gsd.nfreqs for n in calls)


def test_empty_selectors_uses_full_read(long_acq):
    _assert_same(
        read_acq_to_gsdata(long_acq, selectors={}), read_acq_to_gsdata(long_acq)
    )


def test_unknown_selector(long_acq):
    with pytest.raises(ValueError, match="Unrecognized selectors"):
        read_acq_to_gsdata(long_acq, selectors={"bad_selector": {}})


def test_selection_matches_nothing(long_acq):
    with pytest.raises(ACQError, match="Selection matched no data"):
        read_acq_to_gsdata(
            long_acq,
            selectors={
                "time_selector": {"time_range": (T0 - 2 * un.day, T0 - 1 * un.day)}
            },
        )

    with pytest.raises(ACQError, match="matched no frequency channels"):
        read_acq_to_gsdata(
            long_acq,
            selectors={"freq_selector": {"freq_range": (500 * un.MHz, 600 * un.MHz)}},
        )


def test_write_selected_gsh5(long_acq, tmp_path):
    gsd = read_acq_to_gsdata(
        long_acq,
        selectors={
            "freq_selector": {"freq_range": (50 * un.MHz, 100 * un.MHz)},
            "lst_selector": {"lst_range": (6, 12)},
        },
    )
    gsd.write_gsh5(tmp_path / "selected.gsh5")
    new = GSData.from_file(tmp_path / "selected.gsh5")
    np.testing.assert_array_equal(new.data, gsd.data)
