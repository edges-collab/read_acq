"""Tests of writing the data."""
from pathlib import Path

import pytest
from read_acq.gsdata import read_acq_to_gsdata


@pytest.fixture(scope="module")
def sample_acq() -> Path:
    return Path(__file__).parent / "data/sample.acq"


def test_read_as_gsh5(sample_acq: Path):
    gsd = read_acq_to_gsdata(sample_acq, telescope_location="edges")
    assert gsd.filename == sample_acq
