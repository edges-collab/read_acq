#! /usr/bin/env python
"""
Creates an array file for use by sensitivity.py.  The main product is the uv coverage produced by
the array during the
time it takes the sky to drift through the primary beam; other array parameters are also saved.
Array specific information comes from an aipy cal file.  If track is set, produces the uv coverage
for the length specified instead of that set by the primary beam.
"""
from __future__ import division
from __future__ import print_function

import glob

import click

from src import read_acq

main = click.Group()


@main.command()
@click.argument(
    "infile", nargs=-1,
)
@click.option(
    "-o", "--outfile", default=None, type=click.Path(exists=False, dir_okay=False),
)
@click.option(
    "-f",
    "--format",
    default=("mat",),
    multiple=True,
    type=click.Choice(["mat", "h5", "npz"]),
)
def convert(infile, outfile, format):
    fls = []
    for fl in infile:
        fls += glob.glob(fl)
    fls = list(set(fls))

    read_acq.decode_files(fls, outfile=outfile, write_formats=format)