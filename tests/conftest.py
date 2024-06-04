"""Test configuration."""

from pathlib import Path

import pytest


@pytest.fixture(scope="module")
def sample_acq() -> Path:
    return Path(__file__).parent / "data/sample.acq"
