#! /usr/bin/env python
"""
Creates an array file for use by sensitivity.py.  The main product is the uv coverage produced by the array during the
time it takes the sky to drift through the primary beam; other array parameters are also saved.
Array specific information comes from an aipy cal file.  If track is set, produces the uv coverage
for the length specified instead of that set by the primary beam.
"""
from __future__ import division
from __future__ import print_function

import glob

import click

import read_acq

main = click.Group()


@main.command()
@click.argument(
    'infile', nargs=-1,
)
@click.option(
    '-o', '--outfile', default=None,
    type=click.Path(exists=False, dir_okay=False),
)
@click.option(
    '--tcal', default=None,
    help="calibration temperature"
)
@click.option(
    '--tload', default=300,
    help="temperature of load"
)
@click.option(
    '--nchannels', default=16384*2,
    help="number of channels"
)
def convert(infile, outfile, tcal, tload, nchannels):
    fls = []
    for fl in infile:
        fls += glob.glob(fl)

    read_acq.decode_files(fls, outfile=outfile, tcal=tcal, nchannels=nchannels,
                          tload=tload)
