"""Tests of the metadata-only reader, read_metadata."""

from __future__ import annotations

import contextlib
import io
import re
import sys
import warnings
from collections.abc import Callable
from pathlib import Path

import numpy as np
import pytest

import read_acq
from read_acq import decode_file, encode, read_metadata
from read_acq.read_acq import (
    ACQError,
    ACQLineError,
    Ancillary,
    CommentLine,
    DataLine,
    _index_file,
    _iter_entries_without_spectra,
    _read_spectra,
    _Reader,
)

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


def _starts_mid_cycle_with_truncated_leading_spectrum(h, e):
    # decode_file does not read the leading (swpos 1, 2) entries, so this is not an
    # error or a warning.
    e[1][1] = _truncate_spectrum(e[1][1])
    return _join(h, e[1:])


# The first line of the file is read for the metadata, so these modify the second
# leading entry.


def _starts_mid_cycle_with_leading_data_line_missing(h, e):
    return _join(h, [e[1], [e[2][0]], *e[3:]])


def _starts_mid_cycle_with_leading_data_line_cut_before_spectrum(h, e):
    e[2][1] = e[2][1].split(" spectrum ")[0] + "\n"
    return _join(h, e[1:])


def _starts_mid_cycle_with_leading_spectrum_too_long(h, e):
    e[2][1] = e[2][1].rstrip("\n") + "AAAAAAAA\n"
    return _join(h, e[1:])


# Longer than the small read made at each entry, so that the reader has to read more
# to find where the spectrum starts.
_LONG = " " * 600


def _long_comment_line_mid_file(h, e):
    e[7][0] = e[7][0].replace("adcmax ", f"adcmax {_LONG}", 1)
    return _join(h, e)


def _long_front_matter_mid_file(h, e):
    e[7][1] = e[7][1].replace(" spectrum ", f" spectrum {_LONG}", 1)
    return _join(h, e)


def _long_first_head(h, e):
    e[0][0] = e[0][0].replace("adcmax ", f"adcmax {_LONG}", 1)
    e[0][1] = e[0][1].replace(" spectrum ", f" spectrum {_LONG}", 1)
    return _join(h, e)


def _long_leading_head(h, e):
    e[1][0] = e[1][0].replace("adcmax ", f"adcmax {_LONG}", 1)
    e[1][1] = e[1][1].replace(" spectrum ", f" spectrum {_LONG}", 1)
    return _join(h, e[1:])


def _long_header(h, e):
    # Longer than the first read of the file.
    return _join([*h, *(f";--note{i}: {'x' * 60}\n" for i in range(200))], e)


def _header_item_without_value(h, e):
    return _join([*h, ";--novalue\n"], e)


def _junk_line_before_first_comment(h, e):
    return _join([*h, "junk\n"], e)


def _crlf_starting_mid_cycle(h, e):
    return _crlf(h, e[1:])


def _crlf_header_line_split_by_first_read(h, e):
    # The "\r" of a "\r\n" is the last byte of the first read of the file.
    prefix = ";--pad: "
    pad = "x" * (read_acq.read_acq._HEAD_READ_SIZE - len(prefix) - 1)
    return _join([f"{prefix}{pad}\r\n", *h], e)


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
        _starts_mid_cycle_with_truncated_leading_spectrum,
        _starts_mid_cycle_with_leading_data_line_missing,
        _starts_mid_cycle_with_leading_data_line_cut_before_spectrum,
        _starts_mid_cycle_with_leading_spectrum_too_long,
        _long_comment_line_mid_file,
        _long_front_matter_mid_file,
        _long_first_head,
        _long_leading_head,
        _long_header,
        _header_item_without_value,
        _junk_line_before_first_comment,
        _crlf_starting_mid_cycle,
        _crlf_header_line_split_by_first_read,
    ]
}


@pytest.fixture(scope="module")
def clean_file(tmp_path_factory) -> Path:
    return _write_synthetic(tmp_path_factory.mktemp("meta") / "clean.acq")


def _recording_warnings(func, *args):
    with warnings.catch_warnings(record=True) as ws:
        warnings.simplefilter("always")
        out = func(*args)
    return out, [(w.category, str(w.message)) for w in ws]


def _assert_matches_decode_file(path: Path):
    """Check read_metadata and _index_file against decode_file.

    Both must give the same ancillary data and warnings, and the offsets from
    _index_file must point at the spectra that decode_file decodes.
    """
    (_, p, anc), w_decode = _recording_warnings(decode_file, path, False)
    (meta, ancillary), w_meta = _recording_warnings(read_metadata, path)
    (anc_idx, offsets), w_index = _recording_warnings(_index_file, path)

    for m, a in [(meta, ancillary), (anc_idx.meta, anc_idx.data)]:
        assert m == anc.meta
        assert a.keys() == anc.data.keys()
        for key, val in anc.data.items():
            assert a[key].dtype == val.dtype, key
            np.testing.assert_array_equal(a[key], val, err_msg=key)

    assert w_meta == w_decode
    assert w_index == w_decode

    ncycles = len(anc.data["times"])
    assert offsets.shape == (ncycles, 3)
    if ncycles:
        spectra = _read_spectra(path, offsets, anc.meta["nfreq"])
        np.testing.assert_array_equal(spectra, np.transpose(p, (0, 2, 1)))

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


def test_bad_entry_is_yielded_as_none(clean_file: Path, tmp_path: Path):
    header, entries = _split(clean_file)
    path = tmp_path / "bad.acq"
    path.write_bytes(_truncated_spectrum_mid_file(header, entries).encode())

    with (
        path.open("rb", buffering=0) as fl,
        pytest.warns(UserWarning, match="nspec and length of spectrum do not match"),
    ):
        out = list(_iter_entries_without_spectra(path, _Reader(fl), fastspec=True))

    assert len(out) == len(entries)
    for i, item in enumerate(out):
        if i == 7:
            assert item is None
        else:
            entry, offset = item
            assert entry.data.time == entries[i][1][:17]
            assert path.read_bytes()[offset:].startswith(entries[i][1][:17].encode())


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


@pytest.mark.parametrize(
    "reader", [lambda p: decode_file(p, progress=False), read_metadata, _index_file]
)
def test_bad_first_entry_raises(reader, clean_file: Path, tmp_path: Path):
    """As in decode_file, a bad first entry is an error, not just a warning."""
    header, entries = _split(clean_file)
    entries[0][1] = _truncate_spectrum(entries[0][1])
    path = tmp_path / "bad_first.acq"
    path.write_bytes(_join(header, entries).encode())

    with pytest.raises(ACQLineError, match="nspec and length of spectrum"):
        reader(path)


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


# --- A bad first entry is an error, as in decode_file ----------------------------


def _outcome(func, path: Path):
    """Return what ``func(path)`` raises (type and message), and its warnings."""
    with warnings.catch_warnings(record=True) as ws:
        warnings.simplefilter("always")
        try:
            func(path)
        except Exception as e:  # noqa: BLE001 (compare whatever is raised)
            raised = (type(e), str(e))
        else:
            raised = None
    return raised, [(w.category, str(w.message)) for w in ws]


def _first_spectrum_truncated(h, e):
    e[0][1] = _truncate_spectrum(e[0][1])
    return _join(h, e)


def _first_spectrum_too_long(h, e):
    e[0][1] = e[0][1].rstrip("\n") + "AAAAAAAA\n"
    return _join(h, e)


def _first_data_line_cut_before_spectrum(h, e):
    e[0][1] = e[0][1].split(" spectrum ")[0] + "\n"
    return _join(h, e)


def _first_data_line_swpos_mismatch(h, e):
    e[0][1] = f"{e[0][1][:18]}1{e[0][1][19:]}"
    return _join(h, e)


def _first_front_matter_malformed(h, e):
    e[0][1] = "not a time" + e[0][1][17:]
    return _join(h, e)


def _first_comment_malformed(h, e):
    e[0][0] = "# garbage\n"
    return _join(h, e)


def _first_spectrum_truncated_after_leading_entries(h, e):
    e[3][1] = _truncate_spectrum(e[3][1])
    return _join(h, e[1:])


def _leading_comment_malformed(h, e):
    # Not the first comment line, which is read for the metadata.
    e[2][0] = "# garbage\n"
    return _join(h, e[1:])


def _leading_data_line_read_as_comment_line(h, e):
    # decode_file reads this data line as the first swpos=0 comment line, so the
    # comment line after it is read as a data line.
    e[2][1] = f"{e[3][0].rstrip()} spectrum {'A' * 4 * NFREQ}\n"
    return _join(h, e[1:])


def _first_comment_is_last_line(h, e):
    return _join(h, []) + e[0][0]


def _no_comment_lines(h, e):
    return _join(h, [])


def _first_data_line_nul_padded(h, e):
    return _join(h, []) + e[0][0] + _NUL


BAD_FIRST_ENTRY: dict[str, Callable] = {
    f.__name__.lstrip("_"): f
    for f in [
        _first_spectrum_truncated,
        _first_spectrum_too_long,
        _first_data_line_cut_before_spectrum,
        _first_data_line_swpos_mismatch,
        _first_front_matter_malformed,
        _first_comment_malformed,
        _first_spectrum_truncated_after_leading_entries,
        _leading_comment_malformed,
        _leading_data_line_read_as_comment_line,
        _first_comment_is_last_line,
        _no_comment_lines,
        _first_data_line_nul_padded,
    ]
}


@pytest.mark.parametrize("modifier", BAD_FIRST_ENTRY.keys())
def test_bad_first_entry_raises_as_decode_file(
    modifier, clean_file: Path, tmp_path: Path
):
    header, entries = _split(clean_file)
    path = tmp_path / f"{modifier}.acq"
    path.write_bytes(BAD_FIRST_ENTRY[modifier](header, entries).encode())

    expected = _outcome(lambda p: decode_file(p, progress=False), path)
    assert expected[0] is not None, "decode_file should raise"
    assert _outcome(read_metadata, path) == expected
    assert _outcome(_index_file, path) == expected


# --- The header is parsed exactly as it was when read in text mode ---------------


class _TextModeHeader:
    """The text-mode header parsing of Ancillary before #160, as a reference."""

    header_char = ";"

    def __init__(self, fname: Path):
        self.fastspec_version = self._get_fastspec_version(fname)
        self.meta = self.read_metadata(fname)

    def _get_fastspec_version(self, fname: Path):
        with fname.open("r") as fl:
            first_line = fl.readline()
        if first_line.startswith(self.header_char):
            return first_line.split("FASTSPEC")[-1]
        return None

    def _read_header(self, fname: Path):
        out = {}
        name_pattern = re.compile(r"[a-zA-Z_]+")
        with fname.open("r") as fl:
            for line in fl:
                if not line.startswith(self.header_char):
                    break
                if line.startswith("; FASTSPEC"):
                    name = "fastspec_version"
                    val = line.split()[-1]
                else:
                    try:
                        name, val = line.split(": ")
                    except ValueError:
                        warnings.warn(
                            f"In file {fname}, item {line} has no value", stacklevel=1
                        )
                        name = line.split(":")[0]
                        val = ""
                    name = name_pattern.findall(name)[0]
                for tp in [int, float, str]:
                    try:
                        out[name] = tp(val.split()[0])
                        break
                    except IndexError:
                        with contextlib.suppress(ValueError):
                            out[name] = tp(val)
                    except ValueError:
                        pass
        return out

    def read_metadata(self, fname: Path):
        out = self._read_header(fname)
        with fname.open("r") as fl:
            for line in fl:
                if line.startswith("#"):
                    comment = CommentLine.read(line)
                    data = DataLine.read(next(fl), read_spectrum=False)
                    break
            else:
                raise ACQError(f"No comment line found in file {fname}.")
        out.update(
            {
                "temperature": comment.temp,
                "nblk": comment.nblk,
                "nfreq": comment.nspec,
                "freq_min": data.freqmin,
                "freq_max": data.freqmax,
                "freq_res": data.deltaf,
            }
        )
        if comment.resolution is not None:
            out["resolution"] = comment.resolution
        if comment.data_drops is not None:
            out["data_drops"] = comment.data_drops
        return out


def _header_outcome(cls, path: Path):
    out = {}

    def _read(p):
        anc = cls(p)
        out.update(version=anc.fastspec_version, meta=anc.meta)

    return (*_outcome(_read, path), out)


def _fastspec_version_line(h, e):
    return _join(["; FASTSPEC 1.2.3\n", *h], e)


def _fastspec_version_line_only(h, e):
    return "; FASTSPEC 1.2.3\n"


def _empty(h, e):
    return ""


def _non_ascii_header(h, e):
    return _join([*h, *(f";--note{i}: {'µ' * 40}\n" for i in range(100))], e)


def _non_ascii_header_shifted(h, e):
    # Shift the header by a byte, so that one of the two cases splits a multi-byte
    # character at the end of the first read of the file.
    return _join([";\n", *h, *(f";--note{i}: {'µ' * 40}\n" for i in range(100))], e)


def _header_name_without_letters(h, e):
    return _join([*h, ";--123: 4\n"], e)


def _header_value_types(h, e):
    return _join([*h, ";--a: 1\n", ";--b: 2.5 MHz\n", ";--c: text\n", ";--d: \n"], e)


def _cr_only_in_header(h, e):
    return _join([";--a: 1\r;--b: 2\n", *h], e)


HEADER_CASES: dict[str, Callable] = {
    **MODIFIERS,
    **BAD_FIRST_ENTRY,
    **{
        f.__name__.lstrip("_"): f
        for f in [
            _fastspec_version_line,
            _fastspec_version_line_only,
            _empty,
            _non_ascii_header,
            _non_ascii_header_shifted,
            _header_name_without_letters,
            _header_value_types,
            _cr_only_in_header,
        ]
    },
}


@pytest.mark.parametrize("modifier", HEADER_CASES.keys())
def test_header_matches_text_mode_read(modifier, clean_file: Path, tmp_path: Path):
    header, entries = _split(clean_file)
    path = tmp_path / f"{modifier}.acq"
    path.write_bytes(HEADER_CASES[modifier](header, entries).encode())
    assert _header_outcome(Ancillary, path) == _header_outcome(_TextModeHeader, path)


@pytest.mark.parametrize("fname", ["sample.acq", "pxspec.acq"])
def test_header_matches_text_mode_read_real_data(fname):
    path = DATA / fname
    assert _header_outcome(Ancillary, path) == _header_outcome(_TextModeHeader, path)


def test_ancillary_methods_unchanged(clean_file: Path):
    anc = Ancillary(clean_file)
    ref = _TextModeHeader(clean_file)
    assert anc.read_metadata(clean_file) == ref.meta
    assert anc._read_header(clean_file) == ref._read_header(clean_file)
    assert anc._get_fastspec_version(clean_file) == ref.fastspec_version


# --- I/O: how often the file is opened, and how much of it is read ---------------


class _TrickleFile(io.BytesIO):
    """A file whose reads return at most a few bytes, as unbuffered reads may."""

    def read(self, size=-1):
        return super().read(min(size, 7))


@pytest.mark.parametrize(("pos", "size"), [(0, 20), (90, 20), (100, 20), (95, 5)])
def test_reader_reads_all_bytes_asked_for(pos, size):
    data = bytes(range(100))
    reader = _Reader(_TrickleFile(data))
    assert reader.read(pos, size) == (data[pos : pos + size], pos + size > len(data))


_opened: list[str] | None = None


def _audit(event, args):
    if _opened is not None and event == "open" and isinstance(args[0], str | Path):
        _opened.append(str(args[0]))


sys.addaudithook(_audit)


def _count_opens(func, path: Path) -> int:
    global _opened
    _opened = []
    try:
        func(path)
        return sum(Path(p) == path for p in _opened)
    finally:
        _opened = None


@pytest.mark.parametrize(
    "func",
    [
        read_metadata,
        Ancillary,
        lambda p: decode_file(p, progress=False),
        _index_file,
    ],
    ids=["read_metadata", "Ancillary", "decode_file", "_index_file"],
)
@pytest.mark.parametrize("modifier", ["clean", "starts_mid_cycle"])
def test_opens_file_once(func, modifier, clean_file: Path, tmp_path: Path):
    header, entries = _split(clean_file)
    path = tmp_path / f"{modifier}.acq"
    path.write_bytes(MODIFIERS[modifier](header, entries).encode())
    assert _count_opens(func, path) == 1


def _io_counters() -> tuple[int, int]:
    """Return the number of read syscalls and bytes read by this process so far."""
    vals = dict(
        line.split(": ")
        for line in Path("/proc/self/io").read_text().split("\n")
        if line
    )
    return int(vals["syscr"]), int(vals["rchar"])


_needs_proc_io = pytest.mark.skipif(
    not sys.platform.startswith("linux") or not Path("/proc/self/io").exists(),
    reason="needs /proc/self/io to count reads",
)

BIG_NTIMES = 40
BIG_NFREQ = 32768


@pytest.fixture(scope="module")
def big_file(tmp_path_factory) -> Path:
    return _write_synthetic(
        tmp_path_factory.mktemp("big") / "big.acq", ntimes=BIG_NTIMES, nfreq=BIG_NFREQ
    )


# Both read the file by seeking over the spectra, so both are checked for their I/O.
_seeking_readers = pytest.mark.parametrize(
    "reader",
    [lambda p: read_metadata(p)[1], lambda p: _index_file(p)[0].data],
    ids=["read_metadata", "_index_file"],
)


def _measure(reader, path: Path):
    before = _io_counters()
    ancillary = reader(path)
    after = _io_counters()
    return ancillary, after[0] - before[0], after[1] - before[1]


@_needs_proc_io
@_seeking_readers
@pytest.mark.parametrize(
    ("modifier", "ncycles"),
    [("clean", BIG_NTIMES), ("starts_mid_cycle", BIG_NTIMES - 1)],
)
def test_reads_only_entry_heads(
    reader, modifier, ncycles, big_file: Path, tmp_path: Path
):
    """Only a small read per entry: no spectrum is read, not even the first."""
    header, entries = _split(big_file)
    path = tmp_path / f"{modifier}.acq"
    path.write_bytes(MODIFIERS[modifier](header, entries).encode())
    nentries = len(_split(path)[1])

    ancillary, nreads, nbytes = _measure(reader, path)

    assert len(ancillary["times"]) == ncycles
    # A few reads for the header, then one per entry.
    assert nreads <= nentries + 4
    # Each entry's head (comment line and data-line front matter) is ~160 bytes.
    assert nbytes <= 300 * nentries + 8192, f"read {nbytes} bytes"


@_needs_proc_io
@_seeking_readers
def test_long_heads_are_not_read_in_full(reader, big_file: Path, tmp_path: Path):
    """Heads that don't fit in the small read need a larger one, not a full read."""
    header, entries = _split(big_file)
    for entry in entries:
        entry[1] = entry[1].replace(" spectrum ", f" spectrum {_LONG}", 1)
    path = tmp_path / "long_heads.acq"
    path.write_bytes(_join(header, entries).encode())

    ancillary, _, nbytes = _measure(reader, path)

    assert len(ancillary["times"]) == BIG_NTIMES
    assert nbytes < 0.1 * path.stat().st_size, f"read {nbytes} bytes"


def test_public_api():
    assert "read_metadata" in read_acq.__all__
    assert read_acq.read_metadata is read_metadata
