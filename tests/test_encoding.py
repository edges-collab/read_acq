"""Tests of writing ACQ files."""

from pathlib import Path

import numpy as np
import pytest

from read_acq import decode_file, encode
from read_acq.codec import _decode_line, _encode, _encode_line
from read_acq.read_acq import (
    ACQError,
    ACQLineError,
    Ancillary,
    CommentLine,
    DataEntry,
    DataLine,
)


def test_roundtrip():
    a = np.logspace(-5, -1, 20)

    b = _decode_line(_encode(a))
    assert np.allclose(a, b)


def test_roundtrip_ratio():
    a = np.logspace(-5, -1, 20)
    b = np.logspace(-6, -2, 20)

    aa = _decode_line(_encode_line(a, 1000))
    bb = _decode_line(_encode_line(b, 1000))

    assert np.allclose(a[10:] / b[10:], aa[10:] / bb[10:])


@pytest.fixture(scope="module")
def sample_fastspec_acq():
    return Path(__file__).parent / "data/sample.acq"


@pytest.fixture(scope="module")
def sample_pxspec_acq():
    """First six lines from a real file.

    /data5/edges/data/2014_February_Boolardy/mro/low/2016/2016_250_02.acq
    """
    return Path(__file__).parent / "data/pxspec.acq"


@pytest.fixture(scope="module")
def sample_acq_file(tmp_path_factory):
    data = np.logspace(-5, -1, 1500).reshape((3, 5, 100))
    meta = {
        "temperature": 0,
        "nblk": 1000,
        "nfreq": 100,
        "freq_min": 0,
        "freq_max": 200,
        "freq_res": 2,
    }

    fname = tmp_path_factory.mktemp("direc") / "tempfile.acq"

    ancillary = {
        "adcmin": np.ones((10, 3)) * 0.2,
        "adcmax": np.ones((10, 3)) * 0.3,
        "data_drops": np.zeros((10, 3), dtype=int),
        "times": np.array(["2016:080:01:01:01"] * 10),
    }

    encode(fname, data, meta, ancillary)
    return fname


@pytest.fixture(scope="module")
def incomplete_file(sample_acq_file: Path):
    with sample_acq_file.open("r") as fl:
        lines = fl.readlines()
    last_line = lines[-1]
    meta, spec = last_line.split(" spectrum ")
    spec = _decode_line(spec.lstrip())
    spec = spec[: spec.size // 2]
    spec = _encode_line(spec, nblk=1000)
    lines[-1] = f"{meta} spectrum {spec}" + "\n"

    new_file = sample_acq_file.parent / "incomplete.acq"
    with new_file.open("w") as fl:
        fl.writelines(lines)
    return new_file


@pytest.fixture(scope="module")
def incomplete_specline_file(sample_acq_file: Path):
    with sample_acq_file.open("r") as fl:
        lines = fl.readlines()
    last_line = lines[-1]
    meta, spec = last_line.split(" spectrum ")
    spec = _decode_line(spec.lstrip())
    spec = spec[: spec.size // 2]
    spec = _encode_line(spec, nblk=1000)
    lines[-1] = f"{meta} spec"

    new_file = sample_acq_file.parent / "incomplete_spec.acq"
    with new_file.open("w") as fl:
        fl.writelines(lines)
    return new_file


def test_full_file_encode(sample_acq_file: Path):
    with sample_acq_file.open("r") as fl:
        assert fl.readline() == ";--freq_min: 0\n"


def test_incomplete_file(incomplete_file: Path):
    q, p, _ = decode_file(incomplete_file, progress=False)

    # We lost an integration!
    for pp in p:
        assert pp.T.shape == (4, 100)
    assert q.T.shape == (4, 100)


def test_incomplete_specline_file(incomplete_specline_file: Path):
    with incomplete_specline_file.open("r") as fl:
        badline = fl.readlines()[-1]

    with pytest.raises(ACQLineError, match="Could not parse line"):
        DataLine.read(badline, read_spectrum=True)

    with pytest.warns(UserWarning, match="Could not parse line"):
        q, p, _ = decode_file(incomplete_specline_file, progress=False)

    # We lost an integration!
    for pp in p:
        assert pp.T.shape == (4, 100)
    assert q.T.shape == (4, 100)


@pytest.fixture(scope="module")
def pxspec_comment_line() -> str:
    """Return a real comment-line from a pxspec-output file.

    Taken from /data5/edges/data/2014_February_Boolardy/mro/low/2016/2016_250_02.acq
    """
    return (
        "# swpos 0 resolution   24.414 adcmax  0.28123 adcmin -0.31032 temp  0 C "
        "nblk 40960 nspec 32768"
    )


@pytest.fixture(scope="module")
def data_line_prefix() -> str:
    """Return a real data-line prefix from a pxspec-output file.

    Doesn't include the actual spectrum here.
    Taken from /data5/edges/data/2014_February_Boolardy/mro/low/2016/2016_250_02.acq
    """
    return (
        "2016:250:02:06:07 0    0.000 0.006104  200.000  0.3 spectrum /fabricated-data/"
    )


@pytest.fixture(scope="module")
def fastspec_comment_line() -> str:
    """Return a real comment-line from a fastspec-output file.

    From 2014_February_Boolardy/mro/low2/2022/2022_010_16_12_17_low2.acq,
    but manually changed swpos to 1 to enable some "broken" tests below.
    """
    return (
        "# swpos 1 data_drops        0 adcmax  0.21106 adcmin -0.22217 temp  0 C "
        "nblk 65532 nspec 32768"
    )


def test_pxspec_comment_line(pxspec_comment_line):
    c = CommentLine.read(pxspec_comment_line)
    assert c.swpos == 0
    assert c.resolution == 24.414
    assert c.adcmax == 0.28123
    assert c.adcmin == -0.31032
    assert c.temp == 0
    assert c.nblk == 40960
    assert c.nspec == 32768
    assert c.data_drops is None


def test_fastspec_comment_line(fastspec_comment_line):
    c = CommentLine.read(fastspec_comment_line)
    assert c.swpos == 1
    assert c.data_drops == 0
    assert c.adcmax == 0.21106
    assert c.adcmin == -0.22217
    assert c.temp == 0
    assert c.nblk == 65532
    assert c.nspec == 32768
    assert c.resolution is None

    with pytest.raises(ACQError, match="Could not parse line"):
        CommentLine.read(fastspec_comment_line, fastspec=False)


def test_bad_comment_line():
    with pytest.raises(ACQError, match="Could not parse line"):
        CommentLine.read("this is not a comment line")


def test_data_line_prefix_read(data_line_prefix):
    d = DataLine.read(data_line_prefix, read_spectrum=False)
    assert d.time == "2016:250:02:06:07"
    assert d.swpos == 0
    assert d.freqmin == 0
    assert d.freqmax == 200
    assert d.deltaf == 0.006104
    assert d.thing == 0.3


def test_bad_data_line():
    with pytest.raises(ACQError, match="Could not parse line"):
        DataLine.read("this is not a data line spectrum  ", read_spectrum=False)


def test_read_data_entry(pxspec_comment_line, data_line_prefix):
    data_entry = DataEntry.read(
        (pxspec_comment_line, data_line_prefix), read_spectrum=False
    )
    assert data_entry.comment.swpos == 0
    assert data_entry.data.time == "2016:250:02:06:07"


def test_read_data_entry_mismatching_swpos(fastspec_comment_line, data_line_prefix):
    with pytest.raises(ACQLineError, match="swpos of comment and data do not match"):
        DataEntry.read((fastspec_comment_line, data_line_prefix), read_spectrum=False)


def test_read_data_entry_wrong_nspec(sample_fastspec_acq: Path):
    # Get the two data-entry lines from the front of the file.
    lines = []
    with sample_fastspec_acq.open("r") as fl:
        for line in fl:
            if line.startswith("# "):
                lines.append(line)
                lines.append(next(fl))
                break

    lines[0] = lines[0].replace("nspec 32768", "nspec 32767")
    with pytest.raises(ACQLineError, match="nspec and length of spectrum do not match"):
        DataEntry.read(lines)


def test_ancillary_pxspec(sample_pxspec_acq: Path):
    anc = Ancillary(sample_pxspec_acq)
    assert anc.fastspec_version is None
    assert "data_drops" not in anc.meta
    assert anc.meta["resolution"] == 24.414


def test_ancillary_append_bad(sample_pxspec_acq: Path):
    anc = Ancillary(sample_pxspec_acq)

    # Read first two lines of data
    with sample_pxspec_acq.open("r") as fl:
        lines = [fl.readline(), fl.readline()]

    data_entry = DataEntry.read(lines)

    with pytest.raises(ValueError, match="swpos of datas must be"):
        anc.append((data_entry,) * 3)


def test_read_file_starting_with_swpos_1(sample_pxspec_acq: Path, tmp_path: Path):
    with sample_pxspec_acq.open("r") as fl:
        lines = fl.readlines()

    # Put the swpos=1 lines at the start of the file.
    lines = lines[2:4] + lines

    new = tmp_path / "order_test.acq"
    with new.open("w") as fl:
        fl.writelines(lines)

    q, p, _ = decode_file(new)

    # we should have just one complete switch cycle
    assert q.shape == (32768, 1)
    assert p[0].shape == (32768, 1)
    assert p[1].shape == (32768, 1)
    assert p[2].shape == (32768, 1)


def test_read_file_with_incomplete_cycle(sample_pxspec_acq: Path, tmp_path: Path):
    with sample_pxspec_acq.open("r") as fl:
        lines = fl.readlines()

    # Make an imcomplete cycle in the middle.
    lines = lines + lines[:4] + lines

    new = tmp_path / "incomplete_test.acq"
    with new.open("w") as fl:
        fl.writelines(lines)

    q, p, _ = decode_file(new)
    # we should have just two complete switch cycles
    assert q.shape == (32768, 2)
    assert p[0].shape == (32768, 2)
    assert p[1].shape == (32768, 2)
    assert p[2].shape == (32768, 2)


def test_read_file_with_incomplete_cycle_starting_bad(
    sample_pxspec_acq: Path, tmp_path: Path
):
    with sample_pxspec_acq.open("r") as fl:
        lines = fl.readlines()

    # Make an imcomplete cycle in the middle. Now the first bad swpos is not zero.
    lines = lines + lines[2:6] + lines

    new = tmp_path / "incomplete_test2.acq"
    with new.open("w") as fl:
        fl.writelines(lines)

    q, p, _ = decode_file(new)

    # we should have just two complete switch cycles
    assert q.shape == (32768, 2)
    assert p[0].shape == (32768, 2)
    assert p[1].shape == (32768, 2)
    assert p[2].shape == (32768, 2)


def test_roundtrip_on_file(sample_acq: Path, tmp_path: Path):
    q, p, anc = decode_file(sample_acq)
    print(anc.data["data_drops"].shape)
    encode(
        tmp_path / "newfile.acq",
        p=[pp.T for pp in p],
        meta=anc.meta,
        ancillary=anc.data,
    )
    _q, _p, _anc = decode_file(tmp_path / "newfile.acq")

    np.testing.assert_allclose(
        q, _q, rtol=3e-2
    )  # TODO: only this high because of MacOS
    np.testing.assert_allclose(p, _p, rtol=3e-2)
