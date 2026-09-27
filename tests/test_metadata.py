"""Tests of the metadata-only reader, read_metadata."""

from __future__ import annotations

import sys
import warnings
from collections.abc import Callable
from pathlib import Path

import numpy as np
import pytest

import read_acq
from read_acq import decode_file, encode, read_metadata

DATA = Path(__file__).parent / "data"

# Larger than the chunk the fast reader grabs at each entry, so that the reader
# really has to seek over each spectrum.
NFREQ = 1024


def _write_synthetic(
    path: Path,
    ntimes: int = 6,
    nfreq: int = NFREQ,
    adcmin: np.ndarray | None = None,
    data_drops: np.ndarray | None = None,
) -> Path:
    rng = np.random.default_rng(1234)
    p = rng.uniform(1e-5, 1e-1, size=(3, ntimes, nfreq))
    meta = {
        "temperature": 0,
        "nblk": 1000,
        "nfreq": nfreq,
        "freq_min": 0,
        "freq_max": 200,
        "freq_res": 2,
    }
    ancillary = {
        "adcmin": (
            adcmin
            if adcmin is not None
            else rng.uniform(-0.4, -0.2, size=(ntimes, 3)).round(5)
        ),
        "adcmax": rng.uniform(0.2, 0.4, size=(ntimes, 3)).round(5),
        "data_drops": (
            data_drops if data_drops is not None else np.zeros((ntimes, 3), dtype=int)
        ),
        "times": np.array(
            [[f"2016:080:01:{i:02d}:{s:02d}" for s in range(3)] for i in range(ntimes)]
        ),
    }
    encode(path, p, meta, ancillary)

    # Normalise to LF line endings regardless of platform, so that the
    # manipulations below know exactly what bytes are on disk.
    path.write_bytes(path.read_text().encode())
    return path


def _split(path: Path) -> tuple[list[str], list[list[str]]]:
    """Split a file into its header lines and its (comment, data) entries."""
    lines = path.read_bytes().decode().splitlines(keepends=True)
    first = next(i for i, line in enumerate(lines) if line.startswith("#"))
    body = lines[first:]
    return lines[:first], [body[i : i + 2] for i in range(0, len(body), 2)]


def _join(header: list[str], entries: list[list[str]]) -> str:
    return "".join(header) + "".join("".join(e) for e in entries)


def _truncate_spectrum(data_line: str, keep: float = 0.5) -> str:
    front, spec = data_line.rstrip("\n").split(" spectrum ")
    n = int(len(spec) * keep) // 4 * 4
    return f"{front} spectrum {spec[:n]}\n"


# --- Modifications of a clean synthetic file ------------------------------------
# Each takes (header, entries) and returns the full text of the new file.


def _clean(h, e):
    return _join(h, e)


def _truncated_last_spectrum(h, e):
    e[-1][1] = _truncate_spectrum(e[-1][1])
    return _join(h, e)


def _last_line_cut_before_spectrum(h, e):
    e[-1][1] = e[-1][1].split(" spectrum ")[0] + " spec"
    return _join(h, e)


def _last_comment_without_data(h, e):
    return _join(h, e[:-1]) + e[-1][0]


def _no_trailing_newline(h, e):
    return _join(h, e).rstrip("\n")


def _truncated_spectrum_mid_file(h, e):
    # Second entry of the third cycle.
    e[7][1] = _truncate_spectrum(e[7][1])
    return _join(h, e)


def _spectrum_too_long_mid_file(h, e):
    e[7][1] = e[7][1].rstrip("\n") + "AAAAAAAA\n"
    return _join(h, e)


def _spectrum_short_by_less_than_one_channel(h, e):
    # Still decodes to nspec channels (length // 4), so must be kept.
    e[7][1] = e[7][1].rstrip("\n")[:-1] + "\n"
    return _join(h, e)


def _data_line_cut_before_spectrum_mid_file(h, e):
    e[7][1] = e[7][1].split(" spectrum ")[0] + "\n"
    return _join(h, e)


def _missing_entry_mid_file(h, e):
    del e[7]
    return _join(h, e)


def _swpos_mismatch_mid_file(h, e):
    # The data line starts with a 17-character time, then the swpos.
    line = e[7][1]
    assert line[17:20] == " 1 "
    e[7][1] = f"{line[:18]}2{line[19:]}"
    return _join(h, e)


def _starts_mid_cycle(h, e):
    return _join(h, e[1:])


def _starts_with_extra_swpos0(h, e):
    # An extra swpos-0 entry at the start: the first cycle is restarted.
    return _join(h, [list(e[0]), *e])


def _wider_front_matter_mid_file(h, e):
    # The front matter of one data line has a different width to the others.
    e[7][1] = e[7][1].replace(" 0.3 spectrum ", " 0.300 spectrum ", 1)
    return _join(h, e)


def _extra_spaces_before_spectrum(h, e):
    e[7][1] = e[7][1].replace(" spectrum ", " spectrum    ", 1)
    return _join(h, e)


def _tab_before_spectrum(h, e):
    # Stripped when decoding, but not a plain space, so not handled by the fast path.
    e[7][1] = e[7][1].replace(" spectrum ", " spectrum \t", 1)
    return _join(h, e)


def _crlf(h, e):
    return _join(h, e).replace("\n", "\r\n")


def _only_partial_cycle(h, e):
    return _join(h, e[:2])


_NUL = "\x00" * 5000


def _nul_padded_after_complete_cycles(h, e):
    return _join(h, e[:6]) + _NUL


def _nul_padded_mid_cycle(h, e):
    return _join(h, e[:7]) + _NUL


def _nul_padded_data_line(h, e):
    return _join(h, e[:7]) + e[7][0] + _NUL


def _malformed_comment_mid_file(h, e):
    # decode_file then also tries (and fails) to read the data line as a comment.
    e[6][0] = "# garbage\n"
    return _join(h, e)


def _crlf_with_data_line_cut_before_spectrum(h, e):
    # The warning shows the start of the line, including the "\r".
    e[7][1] = e[7][1].split(" spectrum ")[0] + "\n"
    return _crlf(h, e)


MODIFIERS: dict[str, Callable] = {
    f.__name__.lstrip("_"): f
    for f in [
        _clean,
        _truncated_last_spectrum,
        _last_line_cut_before_spectrum,
        _last_comment_without_data,
        _no_trailing_newline,
        _truncated_spectrum_mid_file,
        _spectrum_too_long_mid_file,
        _spectrum_short_by_less_than_one_channel,
        _data_line_cut_before_spectrum_mid_file,
        _missing_entry_mid_file,
        _swpos_mismatch_mid_file,
        _starts_mid_cycle,
        _starts_with_extra_swpos0,
        _wider_front_matter_mid_file,
        _extra_spaces_before_spectrum,
        _tab_before_spectrum,
        _crlf,
        _only_partial_cycle,
        _nul_padded_after_complete_cycles,
        _nul_padded_mid_cycle,
        _nul_padded_data_line,
        _malformed_comment_mid_file,
        _crlf_with_data_line_cut_before_spectrum,
    ]
}


@pytest.fixture(scope="module")
def clean_file(tmp_path_factory) -> Path:
    return _write_synthetic(tmp_path_factory.mktemp("meta") / "clean.acq")


def _decode_and_read(path: Path):
    with warnings.catch_warnings(record=True) as w_decode:
        warnings.simplefilter("always")
        _, _, anc = decode_file(path, progress=False)
    with warnings.catch_warnings(record=True) as w_meta:
        warnings.simplefilter("always")
        meta, ancillary = read_metadata(path)
    return anc, (meta, ancillary), w_decode, w_meta


def _assert_matches_decode_file(path: Path):
    anc, (meta, ancillary), w_decode, w_meta = _decode_and_read(path)

    assert meta == anc.meta
    assert ancillary.keys() == anc.data.keys()
    for key, val in anc.data.items():
        assert ancillary[key].dtype == val.dtype, key
        np.testing.assert_array_equal(ancillary[key], val, err_msg=key)

    def _msgs(ws):
        return [(w.category, str(w.message)) for w in ws]

    assert _msgs(w_meta) == _msgs(w_decode)
    return meta, ancillary


@pytest.mark.parametrize("fname", ["sample.acq", "pxspec.acq"])
def test_matches_decode_file_real_data(fname):
    _assert_matches_decode_file(DATA / fname)


@pytest.mark.parametrize("modifier", MODIFIERS.keys())
def test_matches_decode_file_synthetic(modifier, clean_file: Path, tmp_path: Path):
    header, entries = _split(clean_file)
    path = tmp_path / f"{modifier}.acq"
    path.write_bytes(MODIFIERS[modifier](header, entries).encode())
    _assert_matches_decode_file(path)


def test_expected_cycles(clean_file: Path, tmp_path: Path):
    """Check the number of cycles directly, not just agreement with decode_file."""
    _, ancillary = read_metadata(clean_file)
    assert ancillary["times"].shape == (6, 3)
    assert ancillary["times"][2, 1] == b"2016:080:01:02:01"

    header, entries = _split(clean_file)
    path = tmp_path / "bad.acq"
    path.write_bytes(_truncated_spectrum_mid_file(header, entries).encode())
    with pytest.warns(UserWarning, match="nspec and length of spectrum do not match"):
        _, ancillary = read_metadata(path)

    # The third cycle is lost, the rest are kept.
    assert ancillary["times"].shape == (5, 3)
    assert ancillary["times"][2, 0] == b"2016:080:01:03:00"


def _no_swpos0_entries(h, e):
    return _join(h, [entry for entry in e if not entry[0].startswith("# swpos 0")])


def _swpos0_comment_at_end_only(h, e):
    return _join(h, e[1:3]) + e[3][0]


@pytest.mark.parametrize("modifier", [_no_swpos0_entries, _swpos0_comment_at_end_only])
def test_no_complete_first_entry(modifier, clean_file: Path, tmp_path: Path):
    """With no swpos=0 entry to start from, there are no cycles."""
    header, entries = _split(clean_file)
    path = tmp_path / "no_start.acq"
    path.write_bytes(modifier(header, entries).encode())

    meta, ancillary = _assert_matches_decode_file(path)
    assert meta["nfreq"] == NFREQ
    assert ancillary.keys() == {"adcmax", "adcmin", "times", "data_drops"}
    for val in ancillary.values():
        assert len(val) == 0


def test_varying_comment_line_lengths(tmp_path: Path):
    """Comment lines change length with the sign of adcmin and width of data_drops."""
    ntimes = 6
    adcmin = np.array([[-0.3, 0.3, -0.03]] * ntimes)
    data_drops = np.array([[0, 12345, 7]] * ntimes)
    path = _write_synthetic(
        tmp_path / "varying.acq", ntimes=ntimes, adcmin=adcmin, data_drops=data_drops
    )
    _, ancillary = _assert_matches_decode_file(path)
    np.testing.assert_array_equal(ancillary["data_drops"], data_drops)
    np.testing.assert_allclose(ancillary["adcmin"], adcmin)


def test_does_not_decode_spectra(monkeypatch, clean_file: Path):
    def _fail(*args, **kwargs):
        raise AssertionError("spectrum was decoded")

    monkeypatch.setattr(read_acq.read_acq, "_decode_line", _fail)
    meta, ancillary = read_metadata(clean_file)
    assert meta["nfreq"] == NFREQ
    assert len(ancillary["times"]) == 6


def _bytes_read() -> int:
    with Path("/proc/self/io").open() as fl:
        for line in fl:
            if line.startswith("rchar:"):
                return int(line.split()[1])
    raise RuntimeError("no rchar in /proc/self/io")  # pragma: no cover


@pytest.mark.skipif(
    not sys.platform.startswith("linux") or not Path("/proc/self/io").exists(),
    reason="needs /proc/self/io to count bytes read",
)
def test_reads_small_fraction_of_file(tmp_path: Path):
    path = _write_synthetic(tmp_path / "big.acq", ntimes=40, nfreq=32768)
    size = path.stat().st_size

    before = _bytes_read()
    _, ancillary = read_metadata(path)
    nread = _bytes_read() - before

    assert len(ancillary["times"]) == 40
    assert nread < 0.1 * size, f"read {nread} of {size} bytes"


def test_public_api():
    assert "read_metadata" in read_acq.__all__
    assert read_acq.read_metadata is read_metadata
