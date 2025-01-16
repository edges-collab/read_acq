"""Tests of writing the data."""

import os
from pathlib import Path

import numpy as np
import pytest
from pygsdata.select import select_loads

from read_acq import decode_file
from read_acq.gsdata import fast_lst_setter, read_acq_to_gsdata, write_gsdata_to_acq
from read_acq.read_acq import ACQError


@pytest.mark.parametrize("with_fast_lst_setter", [False, True])
def test_roundtrip_gsh5(sample_acq: Path, tmp_path, with_fast_lst_setter):
    _q, _p, _meta = decode_file(sample_acq)
    gsd = read_acq_to_gsdata(
        sample_acq, lst_setter=fast_lst_setter if with_fast_lst_setter else None
    )
    np.testing.assert_allclose(
        _q.T, (gsd.data[0, 0] - gsd.data[1, 0]) / (gsd.data[2, 0] - gsd.data[1, 0])
    )

    write_gsdata_to_acq(gsd, tmp_path / "new.acq")
    q, p, meta = decode_file(tmp_path / "new.acq", progress=False)

    np.testing.assert_allclose(
        q.T,
        (gsd.data[0, 0] - gsd.data[1, 0]) / (gsd.data[2, 0] - gsd.data[1, 0]),
        rtol=3e-2,  # TODO: this is only this high because of MacOS
    )


def test_multifile_read(sample_acq: Path):
    gsd0 = read_acq_to_gsdata(sample_acq)
    gsd = read_acq_to_gsdata([sample_acq, sample_acq])
    assert gsd.ntimes == gsd0.ntimes * 2


def test_non_existent_file():
    with pytest.raises(FileNotFoundError, match="does not exist"):
        read_acq_to_gsdata("non_existent_file.acq")


def test_non_existent_telescope(sample_acq):
    with pytest.raises(
        ValueError, match="telescope must be a Telescope or name of a KNOWN_TELESCOPE"
    ):
        read_acq_to_gsdata(sample_acq, telescope="non_existent_telescope")


def test_write_single_load(sample_acq: Path, tmp_path):
    gsd = select_loads(read_acq_to_gsdata(sample_acq), loads=("ant",))

    with pytest.raises(ValueError, match="Can only encode 3-load data to ACQ file"):
        write_gsdata_to_acq(gsd, tmp_path / "new.acq")


def test_no_files_given():
    with pytest.raises(ValueError, match="No files given to read"):
        read_acq_to_gsdata(path=[])


def test_only_empty_files_given(sample_acq, tmp_path):
    outpath = tmp_path / "empty.acq"

    os.system(f"head -n 31 {sample_acq} > {outpath}")

    gsd = read_acq_to_gsdata([sample_acq])

    with pytest.raises(ACQError, match="No data in any files"):
        read_acq_to_gsdata([outpath])

    with pytest.raises(ACQError, match="No data in any files"):
        read_acq_to_gsdata([outpath, outpath])

    gsd1 = read_acq_to_gsdata([sample_acq, outpath])
    assert gsd1.ntimes == gsd.ntimes

    gsd1 = read_acq_to_gsdata([outpath, sample_acq])
    assert gsd1.ntimes == gsd.ntimes

    gsd1 = read_acq_to_gsdata([outpath, outpath, sample_acq])
    assert gsd1.ntimes == gsd.ntimes

    gsd1 = read_acq_to_gsdata([outpath, sample_acq, outpath])
    assert gsd1.ntimes == gsd.ntimes

    gsd1 = read_acq_to_gsdata([sample_acq, outpath, sample_acq])
    assert gsd1.ntimes == gsd.ntimes * 2
