[![bamAlignCleaner-CI](https://github.com/maxibor/bamAlignCleaner/actions/workflows/ci.yml/badge.svg)](https://github.com/maxibor/bamAlignCleaner/actions/workflows/ci.yml) [![PyPI](https://img.shields.io/pypi/v/bamAlignCleaner?color=green&label=Install%20with%20pip)](https://pypi.org/project/bamAlignCleaner/) [![](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](https://anaconda.org/bioconda/bamAlignCleaner) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6344904.svg)](https://doi.org/10.5281/zenodo.6344904)



# bamAlignCleaner

A simple utility tool to remove references with no aligned reads in a `bam`/`cram` file

## Installation

* with pip 

```bash
pip install bamAlignCleaner 
```

* with conda

```bash
conda install -c bioconda bamAlignCleaner
```

## Usage

* Quick start

```bash
bamAlignCleaner input.bam
```

* Advanced

```bash
bamAlignCleaner --output output.bam --method parse input.bam
```

## Help

```bash
$ bamAlignCleaner --help
Usage: bamAlignCleaner [OPTIONS] BAM

  bamAlignCleaner: removes unaligned references in BAM/CRAM alignment files
  * Homepage: https://github.com/maxibor/bamAlignCleaner
  * Author: Maxime Borry

  BAM: BAM alignment file (sorted, and optionally indexed)

Options:
  --version                       Show the version and exit.
  -m, --method [parse|index_stat]
                                  unaligned reference removal method
                                  [default: parse]
  -o, --output FILE               filtered bam file [default: STDOUT]
  --help                          Show this message and exit.
```

## Methods

bamAlignCleaner uses either one of the two following methods to remove references not having reads mapped against them.

* The `parse` method goes through each read of the alignment file and keeps the references if the read maps to it. *This method should be faster if you have more references than reads.*
* The `check_index` checks index for the number of mapped reads to each reference. *This method should be faster if you have more reads than references.*

## Cite

```bibtex
@software{maxime_borry_2022_6344905,
  author       = {Maxime Borry},
  title        = {maxibor/bamAlignCleaner: bamAlignCleaner},
  month        = mar,
  year         = 2022,
  publisher    = {Zenodo},
  version      = {0.2.2},
  doi          = {10.5281/zenodo.6344904},
  url          = {https://doi.org/10.5281/zenodo.6344904}
}
```