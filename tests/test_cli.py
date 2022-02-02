import h5py
import numpy as np
from click.testing import CliRunner
from pathlib import Path

from read_acq import decode_file
from read_acq.cli import convert


def test_convert(tmp_path_factory):
    cli = CliRunner()

    data = Path(__file__).parent / "data/sample.acq"
    outfile = tmp_path_factory.mktemp("direc") / "tempfile.h5"

    result = cli.invoke(convert, [str(data), "--outfile", str(outfile), "--fmt", "h5"])
    assert result.exit_code == 0
    print(result.output)

    Q, p, meta = decode_file(data, meta=True)

    with h5py.File(outfile, "r") as fl:
        qq = fl["spectra"]["Q"][...]

    assert np.allclose(qq[~np.isnan(qq)], Q[~np.isnan(Q)])
