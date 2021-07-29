[![bamAlignCleaner-CI](https://github.com/maxibor/bamAlignCleaner/actions/workflows/ci.yml/badge.svg)](https://github.com/maxibor/bamAlignCleaner/actions/workflows/ci.yml) [![PyPI](https://img.shields.io/pypi/v/bamAlignCleaner?color=green&label=Install%20with%20pip)](https://pypi.org/project/bamAlignCleaner/) [![Anaconda-Server Badge](https://anaconda.org/maxibor/bamaligncleaner/badges/installer/conda.svg)](https://anaconda.org/maxibor/bamaligncleaner)

# bamAlignCleaner

A simple utility tool to remove references with no aligned reads in a `bam`/`cram` file

## Installation

With pip 
```bash
pip install bamAlignCleaner 
```

With conda

```bash
conda install -c bioconda -c conda-forge -c maxibor bamAlignCleaner
```

## Usage

```bash
bamAlignCleaner -o output.bam input.bam
```
