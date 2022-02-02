#! /usr/bin/env python
"""Command-Line Interface for read_acq."""

import click
import glob
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
    "--fmt",
    default="h5",
    type=click.Choice(["h5", "mat", "npz"]),
)
def convert(infile, outfile, fmt):
    """Convert an ACQ file to a different format."""
    fls = []
    for fl in infile:
        fls += glob.glob(fl)
    fls = list(set(fls))

    for fl in tqdm.tqdm(
        fls, disable=len(fls) < 5, desc="Processing files", unit="files"
    ):
        convert_file(fl, outfile=outfile, write_format=fmt)
