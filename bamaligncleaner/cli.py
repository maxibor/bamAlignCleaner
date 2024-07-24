#!/usr/bin/env python3

import logging
import sys

import click

from bamaligncleaner import __version__
from bamaligncleaner.main import filter_bam


@click.command()
@click.version_option(__version__)
@click.argument("bam", type=click.Path(exists=True))
@click.option(
    "-m",
    "--method",
    type=click.Choice(["parse", "index_stat"]),
    default="parse",
    show_default=True,
    help="unaligned reference removal method",
)
@click.option(
    "-o",
    "--output",
    type=click.Path(writable=True, dir_okay=False, file_okay=True),
    default="-",
    help="filtered bam file [default: STDOUT]",
)
@click.option(
    "-s",
    "--splits",
    default=1,
    show_default=True,
    help="the number of bam files to split the bam file to"
)
def cli(bam, method, output, splits):
    """\b
    bamAlignCleaner: removes unaligned references in BAM/CRAM alignment files
    * Homepage: https://github.com/maxibor/bamAlignCleaner
    * Author: Maxime Borry, Alex Hübner

    BAM: BAM alignment file (sorted, and optionally indexed)
    """
    if splits > 1 and output == "-":
        logging.error("Splitting the BAM file does not work with writing the "
                      "filtered data to STDOUT. Please specify an output "
                      "prefix.")
        sys.exit(1)

    filter_bam(bam, method, output, splits)


if __name__ == "__main__":
    cli()
