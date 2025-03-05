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
    "-r",
    "--reflist",
    type=click.Path(exists=True),
    help="File listing references to keep in output bam",
)
@click.option(
    "-o",
    "--output",
    type=click.Path(writable=True, dir_okay=False, file_okay=True),
    default="-",
    help="filtered bam file [default: STDOUT]. When option --splits > 1, then its the prefix for naming the filtered bam files",
)
@click.option(
    "-s",
    "--splits",
    default=1,
    show_default=True,
    help="the number of bam files to split the bam file to"
)
<<<<<<< HEAD
@click.option(
    "--splitmode",
    type=click.Choice(["contigs", "reads"]),
    default="contigs",
    show_default=True,
    help="method to split the bam file into multiple output files"
)
def cli(bam, method, reflist, output, splits, splitmode):
    """\b
    bamAlignCleaner: removes unaligned references in BAM/CRAM alignment files
    * Homepage: https://github.com/maxibor/bamAlignCleaner
    * Author: Maxime Borry, Alex Hübner

    BAM: BAM alignment file (sorted, and optionally indexed)
    """
    if splits > 100:
        logging.error("It is not supported to split the BAM file into more "
                      "than 100 files.")
        sys.exit(1)
    if splits > 1 and output == "-":
        logging.error("Splitting the BAM file does not work with writing the "
                      "filtered data to STDOUT. Please specify an output "
                      "prefix.")
        sys.exit(1)
    if splits > 1 and method == "index_stat":
        logging.error("Splitting the BAM file does not work with using the "
                      "method index_stat. Please specify the method 'parse'.")
        sys.exit(1)

    filter_bam(bam, method, reflist, output, splits, splitmode)


if __name__ == "__main__":
    cli()
