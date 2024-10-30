"""Tests of the CLI."""

from pathlib import Path

import numpy as np
from click.testing import CliRunner
from pygsdata import GSData

from read_acq import decode_file
from read_acq.cli import convert


def test_convert(tmp_path_factory):
    cli = CliRunner()

    data = Path(__file__).parent / "data/sample.acq"
    outfile = tmp_path_factory.mktemp("direc") / "tempfile.gsh5"
    result = cli.invoke(convert, [str(data), "--outfile", str(outfile)])
    assert result.exit_code == 0

    q, p, meta = decode_file(data)
    gsd = GSData.from_file(outfile)
    assert np.allclose(gsd.data[0, 0].T, p[0])
