"""Benchmark the I/O done by read_metadata (and the header read shared with decode_file).

Writes synthetic ACQ files and, for each reader, reports the wall time, the number of
times the file is opened, and (on Linux) the number of read syscalls and bytes read.
On a network filesystem, per-file cost is dominated by latency, so the opens and
reads matter more than the wall time on local disk.

Usage::

    python benchmarks/bench_metadata.py [--ntimes 40] [--nfreq 32768] [--repeat 5]
"""

from __future__ import annotations

import argparse
import os
import sys
import tempfile
import time
from pathlib import Path

import numpy as np
from bench_read import write_file

from read_acq import decode_file, read_metadata
from read_acq.read_acq import Ancillary, _index_file

_PROC_IO = Path("/proc/self/io")
_opened: list[str] = []


def _audit(event, args):
    if event == "open" and isinstance(args[0], str | bytes | os.PathLike):
        _opened.append(os.fsdecode(args[0]))


def _io_counters() -> tuple[int, int] | None:
    """Return the (read syscalls, bytes read) of this process so far, if available."""
    if not _PROC_IO.exists():
        return None
    vals = dict(line.split(": ") for line in _PROC_IO.read_text().splitlines())
    return int(vals["syscr"]), int(vals["rchar"])


def measure(func, fname: Path, repeat: int) -> dict:
    """Measure one call of ``func(fname)`` (and the best time of ``repeat`` calls)."""
    best = np.inf
    for _ in range(repeat):
        t0 = time.perf_counter()
        func(fname)
        best = min(best, time.perf_counter() - t0)

    _opened.clear()
    before = _io_counters()
    func(fname)
    after = _io_counters()

    out = {"time": best, "opens": sum(p == str(fname) for p in _opened)}
    if before is not None:
        out["reads"] = after[0] - before[0]
        out["bytes"] = after[1] - before[1]
    return out


def _starts_mid_cycle(fname: Path) -> Path:
    """Copy a file, dropping its first entry, so it starts at swpos 1."""
    lines = fname.read_bytes().splitlines(keepends=True)
    first = next(i for i, line in enumerate(lines) if line.startswith(b"#"))
    out = fname.with_name(f"mid_{fname.name}")
    out.write_bytes(b"".join(lines[:first] + lines[first + 2 :]))
    return out


def main():
    """Run the benchmarks."""
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--ntimes", type=int, default=40)
    parser.add_argument("--nfreq", type=int, default=32768)
    parser.add_argument("--repeat", type=int, default=5)
    args = parser.parse_args()

    sys.addaudithook(_audit)
    readers = {
        "read_metadata": read_metadata,
        "Ancillary": Ancillary,
        "_index_file": _index_file,
        "decode_file": lambda f: decode_file(f, progress=False),
    }

    with tempfile.TemporaryDirectory() as tmp:
        clean = Path(tmp) / "bench.acq"
        write_file(clean, args.ntimes, args.nfreq)
        files = {"clean": clean, "starts mid-cycle": _starts_mid_cycle(clean)}

        print(
            f"Files: {args.ntimes} cycles x {args.nfreq} channels, "
            f"{clean.stat().st_size / 1024**2:.1f} MB, {3 * args.ntimes} entries\n"
        )
        print(
            f"{'file':<18}{'reader':<15}{'time (ms)':>10}{'opens':>7}"
            f"{'reads':>8}{'KB read':>10}"
        )
        for label, fname in files.items():
            for name, func in readers.items():
                r = measure(func, fname, args.repeat)
                reads = r.get("reads", "n/a")
                kb = f"{r['bytes'] / 1024:.1f}" if "bytes" in r else "n/a"
                print(
                    f"{label:<18}{name:<15}{1000 * r['time']:>10.2f}{r['opens']:>7}"
                    f"{reads:>8}{kb:>10}"
                )


if __name__ == "__main__":
    main()
