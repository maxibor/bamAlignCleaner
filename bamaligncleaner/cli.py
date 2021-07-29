#!/usr/bin/env python3

import click
from bamaligncleaner import __version__
from bamaligncleaner.main import filter_bam

@click.command()
@click.version_option(__version__)
@click.argument("bam", type=click.Path(exists=True))
@click.option(
    "-o",
    "--output",
    type=click.Path(writable=True, dir_okay=False, file_okay=True),
    default="no_unaligned_ref.bam",
    show_default=True,
    help="filtered bam file",
)
def cli(bam, output):
    """
    BAM: BAM alignment file (indexed and sorted)
    """
    filter_bam(bam, output)


if __name__ == '__main__':
    cli()
