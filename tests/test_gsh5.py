"""Tests of writing the data."""

from pathlib import Path

import numpy as np
from read_acq import decode_file
from read_acq.gsdata import read_acq_to_gsdata, write_gsdata_to_acq


def test_roundtrip_gsh5(sample_acq: Path, tmp_path):
    _q, _p, _meta = decode_file(sample_acq, meta=True)
    gsd = read_acq_to_gsdata(sample_acq, telescope_location="edges")
    np.testing.assert_allclose(
        _q.T, (gsd.data[0, 0] - gsd.data[1, 0]) / (gsd.data[2, 0] - gsd.data[1, 0])
    )

    write_gsdata_to_acq(gsd, tmp_path / "new.acq")
    q, p, meta = decode_file(tmp_path / "new.acq", progress=False, meta=True)

    np.testing.assert_allclose(
        q.T,
        (gsd.data[0, 0] - gsd.data[1, 0]) / (gsd.data[2, 0] - gsd.data[1, 0]),
        rtol=3e-2,  # TODO: this is only this high because of MacOS
    )
