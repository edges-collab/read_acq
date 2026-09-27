"""Check that the examples in the README run."""

import re
import warnings
from pathlib import Path

import numpy as np
import pytest
from astropy import units as un
from astropy.time import Time

from read_acq import decode_file, encode

README = Path(__file__).parent.parent / "README.md"
NTIMES = 24
NFREQ = 512


def _write(path: Path, start: str):
    """Write a synthetic ACQ file spanning a day, so that LST selections match."""
    rng = np.random.default_rng(0)
    t0 = Time(start, format="yday", scale="utc")
    times = t0 + np.arange(NTIMES)[:, None] * un.hour + [0, 13, 26] * un.s
    encode(
        path,
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
            "adcmax": rng.uniform(0, 1, size=(NTIMES, 3)).round(5),
            "adcmin": rng.uniform(-1, 0, size=(NTIMES, 3)).round(5),
        },
    )


def _python_blocks() -> list[str]:
    return re.findall(r"```python\n(.*?)```", README.read_text(), flags=re.DOTALL)


def test_readme_has_examples():
    assert len(_python_blocks()) >= 5


def test_readme_examples_run(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    monkeypatch.chdir(tmp_path)
    _write(tmp_path / "my_data.acq", "2023:070:00:00:00")
    _write(tmp_path / "my_data_2.acq", "2023:071:00:00:00")

    namespace: dict = {}
    with warnings.catch_warnings():
        # e.g. astropy's warnings about IERS tables, which have nothing to do with
        # the examples.
        warnings.simplefilter("ignore")
        for block in _python_blocks():
            exec(block, namespace)  # noqa: S102

    # Check the examples do what the README says they do.
    assert namespace["ncycles"] == NTIMES
    assert namespace["q"].shape == (NFREQ, NTIMES)
    assert 0 < namespace["data"].ntimes < NTIMES

    _, p, _ = decode_file(tmp_path / "copy.acq", progress=False)
    # The encoding stores ~6 significant figures.
    np.testing.assert_allclose(
        p, [namespace["p0"], namespace["p1"], namespace["p2"]], rtol=1e-5
    )
