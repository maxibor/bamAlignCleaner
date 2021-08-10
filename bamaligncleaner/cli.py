#!/usr/bin/env python3

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
def cli(bam, method, output):
    """\b
    bamAlignCleaner: removes unaligned references in BAM/CRAM alignment files
    * Homepage: https://github.com/maxibor/bamAlignCleaner
    * Author: Maxime Borry

    BAM: BAM alignment file (sorted, and optionally indexed)
    """
    filter_bam(bam, method, output)


if __name__ == "__main__":
    cli()
