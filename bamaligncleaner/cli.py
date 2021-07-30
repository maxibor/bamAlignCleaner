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
    type=click.Choice(['parse', 'index_stat']),
    default="parse",
    show_default=True,
    help="removal method. Try slow if the fast method isn't working",
)
@click.option(
    "-o",
    "--output",
    type=click.Path(writable=True, dir_okay=False, file_okay=True),
    default="no_unaligned_ref.bam",
    show_default=True,
    help="filtered bam file",
)
def cli(bam, method, output):
    """
    bamAlignCleaner: removes unaligned references in BAM alignment file

    BAM: BAM alignment file (indexed and sorted)
    """
    filter_bam(bam, method, output)


if __name__ == '__main__':
    cli()
