"""Benchmark full vs. partial reads of an ACQ file.

Writes a synthetic ACQ file (by default 200 switch cycles of 32768 channels, ~80 MB,
spanning ~24 hours) and times reading it in full and with various selectors.

Usage::

    python benchmarks/bench_read.py [--ntimes 200] [--nfreq 32768] [--repeat 3]
"""

from __future__ import annotations

import argparse
import tempfile
import time
import tracemalloc
from pathlib import Path

import numpy as np
from astropy import units as un
from astropy.time import Time

from read_acq import encode, read_acq_to_gsdata


def write_file(fname: Path, ntimes: int, nfreq: int):
    """Write a synthetic ACQ file."""
    rng = np.random.default_rng(0)
    start = Time("2023:070:00:00:00", format="yday", scale="utc")
    step = 24 * un.hour / ntimes
    times = start + np.arange(ntimes)[:, None] * step + [0, 13, 26] * un.s

    encode(
        fname,
        p=list(rng.uniform(1e-3, 1, size=(3, ntimes, nfreq))),
        meta={
            "temperature": 25,
            "nblk": 2974,
            "nfreq": nfreq,
            "freq_min": 0.0,
            "freq_res": 200.0 / nfreq,
            "freq_max": 200.0,
        },
        ancillary={
            "times": np.array(
                [[t.strftime("%Y:%j:%H:%M:%S") for t in row] for row in times]
            ),
            "adcmax": rng.uniform(0, 1, size=(ntimes, 3)),
            "adcmin": rng.uniform(-1, 0, size=(ntimes, 3)),
        },
    )
    return start


def bench(fname: Path, selectors: dict | None, repeat: int) -> tuple[float, float]:
    """Return the best wall time (s) and peak traced memory (MB) of a read."""
    best = np.inf
    for _ in range(repeat):
        t0 = time.perf_counter()
        read_acq_to_gsdata(fname, selectors=selectors)
        best = min(best, time.perf_counter() - t0)

    tracemalloc.start()
    read_acq_to_gsdata(fname, selectors=selectors)
    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    return best, peak / 1024**2


def main():
    """Run the benchmarks."""
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--ntimes", type=int, default=200)
    parser.add_argument("--nfreq", type=int, default=32768)
    parser.add_argument("--repeat", type=int, default=3)
    args = parser.parse_args()

    with tempfile.TemporaryDirectory() as tmp:
        fname = Path(tmp) / "bench.acq"
        t0 = write_file(fname, args.ntimes, args.nfreq)
        print(
            f"File: {args.ntimes} cycles x {args.nfreq} channels, "
            f"{fname.stat().st_size / 1024**2:.1f} MB\n"
        )

        cases = {
            "full read": None,
            "time: first 10%": {
                "time_selector": {"time_range": (t0, t0 + 2.4 * un.hour)}
            },
            "lst: 6h window": {"lst_selector": {"lst_range": (6, 12)}},
            "freq: 50-100 MHz": {
                "freq_selector": {"freq_range": (50 * un.MHz, 100 * un.MHz)}
            },
            "time 10% + freq 50-100": {
                "time_selector": {"time_range": (t0, t0 + 2.4 * un.hour)},
                "freq_selector": {"freq_range": (50 * un.MHz, 100 * un.MHz)},
            },
        }

        print(f"{'case':<26}{'time (s)':>10}{'speedup':>10}{'peak MB':>10}")
        ref = None
        for label, selectors in cases.items():
            t, mem = bench(fname, selectors, args.repeat)
            ref = ref or t
            print(f"{label:<26}{t:>10.3f}{ref / t:>9.1f}x{mem:>10.1f}")


if __name__ == "__main__":
    main()
