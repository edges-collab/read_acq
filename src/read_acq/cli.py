#! /usr/bin/env python
"""
Creates an array file for use by sensitivity.py.  The main product is the uv coverage produced by
the array during the
time it takes the sky to drift through the primary beam; other array parameters are also saved.
Array specific information comes from an aipy cal file.  If track is set, produces the uv coverage
for the length specified instead of that set by the primary beam.
"""

import glob

import click
import tqdm

from . import convert_file

main = click.Group()


@main.command()
@click.argument(
    "infile",
    nargs=-1,
)
@click.option(
    "-o",
    "--outfile",
    default=None,
    type=click.Path(exists=False, dir_okay=False),
)
@click.option(
    "-f",
    "--format",
    default="h5",
    type=click.Choice(["h5", "mat", "npz"]),
)
def convert(infile, outfile, format):
    fls = []
    for fl in infile:
        fls += glob.glob(fl)
    fls = list(set(fls))

    for fl in tqdm.tqdm(
        fls, disable=len(fls) < 5, desc="Processing files", unit="files"
    ):
        convert_file(fl, outfile=outfile, write_format=format)
