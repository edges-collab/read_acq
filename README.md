# read-acq

**Read EDGES ACQ spectrum files.**

[![image](https://travis-ci.org/edges-collab/read_acq.svg?branch=master)](https://travis-ci.org/edges-collab/read_acq)

[![image](https://codecov.io/gh/edges-collab/read_acq/branch/master/graph/badge.svg)](https://travis-ci.org/edges-collabcodecov.io/gh/edges-collab/read_acq)

[![image](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

## Installation

In a new/existing python environment, run `pip install
git+https://github.com/edges-collab/read_acq`.

If you wish to develop `read_acq`, do the following:

    git clone https://github.com/edges-collab/read_acq
    cd read_acq
    pip install -e .

## Usage

`read_acq` can be used from Python or from the command line.

### Which reader do I need?

| You want | Use | Reads from disk |
|---|---|---|
| A [`GSData`](https://github.com/edges-collab/pygsdata) object (what most analysis code works with) | `read_acq_to_gsdata`, or `GSData.from_file` | the spectra you select |
| Only part of a file: a time, LST or frequency range, or one load | `read_acq_to_gsdata(..., selectors=...)` | only the selected spectra and channels |
| Only the metadata and per-cycle ancillary data (times, ADC max/min, data drops) | `read_metadata` | a few hundred bytes per spectrum |
| The raw decoded arrays | `decode_file` | the whole file |

All of these handle a file in the same way (see [Damaged and incomplete
files](#damaged-and-incomplete-files)), and give the same data for it: the cheaper
readers skip over what they don't need rather than reading it.

### Reading into a GSData object

```python
from read_acq import read_acq_to_gsdata

data = read_acq_to_gsdata("my_data.acq", telescope="edges-low")
```

Several files can be given, and they are joined along the time axis:

```python
data = read_acq_to_gsdata(["my_data.acq", "my_data_2.acq"], telescope="edges-low")
```

`read_acq` also registers itself as the reader of `.acq` files for `pygsdata`, so
`GSData.from_file("my_data.acq")` works too. See `help(read_acq_to_gsdata)` for all
the options, e.g. `lst_setter=read_acq.gsdata.fast_lst_setter` for a much faster
approximate calculation of the LSTs (within ~1 s of astropy's for the test data).

### Reading only part of a file

Pass `selectors` to select data while reading. Only the selected spectra, and only
the selected channels of each, are decoded, so this is much faster (and uses much
less memory) than reading everything and then selecting:

```python
from astropy import units as un

data = read_acq_to_gsdata(
    "my_data.acq",
    telescope="edges-low",
    selectors={
        "lst_selector": {"lst_range": (6, 12)},
        "freq_selector": {"freq_range": (50 * un.MHz, 100 * un.MHz)},
    },
)
```

The keys are `freq_selector`, `time_selector`, `lst_selector` and `load_selector`,
and each maps to the keyword arguments of the matching function in
[`pygsdata.select`](https://github.com/edges-collab/pygsdata) (`select_freqs`,
`select_times`, `select_lsts` and `select_loads`). The result is exactly what
reading everything and then applying those functions (in that order) gives.
`GSData.from_file("my_data.acq", selectors=...)` does the same.

### Reading only the metadata

To get only the metadata and the per-cycle ancillary data, use `read_metadata`. It
never reads the spectra: it opens the file once and seeks from one spectrum to the
next, reading only a few hundred bytes around each. This makes it cheap even on a
network filesystem, where each read costs a round trip, so it is a good way to scan
many files (e.g. to find which ones cover a given time range):

```python
from read_acq import read_metadata

meta, ancillary = read_metadata("my_data.acq")
ncycles = len(ancillary["times"])
```

`meta` is a dict of the file's header and settings (`nfreq`, `freq_min`,
`freq_max`, `freq_res`, `nblk`, `temperature`, ...). `ancillary` has an array for
each of `times` (as `YYYY:DDD:HH:MM:SS` byte strings), `adcmax`, `adcmin` and, for
fastspec files, `data_drops`, each with shape `(ncycles, 3)`: one column for each
switch position (antenna, load, load + noise source).

### Reading the raw arrays

`decode_file` decodes the whole file into `numpy` arrays:

```python
from read_acq import decode_file

q, (p0, p1, p2), anc = decode_file("my_data.acq", progress=False)
```

`p0`, `p1` and `p2` are the powers of the three switch positions (antenna, load,
load + noise source), each with shape `(nfreq, ncycles)`, and `q = (p0 - p1) / (p2 -
p1)`. `anc` holds the metadata (`anc.meta`), the per-cycle ancillary data
(`anc.data`, as from `read_metadata`) and the frequencies in MHz
(`anc.frequencies`).

### Damaged and incomplete files

ACQ files are often cut short, e.g. when an acquisition is interrupted. All the
readers handle this in the same way:

- Only complete switch cycles (swpos 0, 1 and 2, in order) are returned. Reading
  starts at the first swpos-0 entry, and a cycle with a missing or bad entry is
  dropped.
- A bad entry after the first (e.g. a spectrum of the wrong length) gives a
  warning, and reading carries on.
- If the rest of a file is NUL bytes (an interrupted write), you get a warning and
  the complete cycles read before it.
- A bad *first* swpos-0 entry, or a malformed comment line before it, is an error,
  as is a file with no comment lines at all.

### Writing ACQ files

`encode(filename, [p0, p1, p2], meta, ancillary)` writes an ACQ file from the powers
of the three switch positions, each with shape `(ncycles, nfreq)` (i.e. transposed
from what `decode_file` returns), and metadata and ancillary data like those from
`read_metadata`. For example, to write a copy of a file (the powers are stored to
about 1 part in 10<sup>6</sup>, so they come back very slightly changed):

```python
from read_acq import encode

encode("copy.acq", [p0.T, p1.T, p2.T], anc.meta, anc.data)
```

`write_gsdata_to_acq` writes a `GSData` object to an ACQ file.

### CLI

The `acq convert` command converts ACQ files to a single GSH5 file (`pygsdata`'s
HDF5 format):

    acq convert data/2023_070*.acq

The input files can include glob-style patterns. Matching files are sorted and
joined along the time axis. By default the output is named by the year and day of
the first time in the data (e.g. `2023_070.gsh5`), and is written in the current
directory. See `acq convert --help` for the options.

## Performance

`benchmarks/` has scripts to measure the readers on synthetic files:

- `python benchmarks/bench_read.py` times full and partial reads with
  `read_acq_to_gsdata`.
- `python benchmarks/bench_metadata.py` reports, for each reader, the wall time, the
  number of times the file is opened and (on Linux) the number of reads and bytes
  read. On a network filesystem these counts, not the wall time on local disk,
  are what decide how long reading takes.
