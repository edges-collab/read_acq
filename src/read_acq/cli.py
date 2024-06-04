#! /usr/bin/env python
"""Command-Line Interface for read_acq."""

from __future__ import annotations

from pathlib import Path

import click

from .gsdata import read_acq_to_gsdata

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
    "-d",
    "--direc",
    default=".",
    type=click.Path(exists=True, file_okay=False),
)
@click.option("--telescope", "--loc", default="edges-low", type=str)
@click.option("-n", "--name", default="{year}_{day}", type=str)
def convert(infile, outfile: str | None, direc: str, telescope: str, name: str):
    """Convert an ACQ file to a .gsh5 format file.

    Multiple input files can be provided, and they will be joined together in a single
    output file. Furthermore, the input files can have glob-style patterns in them so
    that they expand to a list of files, e.g.:

        $ acq convert 2023_070*.acq

    If there is more than one file matching the pattern, the files will be sorted, then
    all read and output as a single GSH5 file (by default, the file "2023_070.gsh5").
    """
    infiles = [Path(fl) for fl in infile]

    fls = []
    for fl in infiles:
        fls += sorted(fl.parent.glob(fl.name))

    gsd = read_acq_to_gsdata(fls, telescope=telescope, name=name)
    if outfile is None:
        outfile = f"{gsd.name}.gsh5"

    gsd.write_gsh5(Path(direc) / outfile)
