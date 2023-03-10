[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](./LICENSE)
[![Maintainer](https://img.shields.io/badge/Maintainer-davidmeijer-blue)](https://github.com/davidmeijer)
[![Generic badge](https://img.shields.io/badge/Version-alpha-green.svg)](https://shields.io/)

# COMET: COmparative Molecular Enrichment Tool

<img src="https://github.com/moltools/COMET/blob/main/logo.png" alt="logo" width="100">

COMET (COmparative Molecular Enrichment Tool) is a Python library and command line tool for extracting comparative enriched substructures from groups of molecules.

## Installing

Clone the repository and move to the root of the repository. Install COMET by running `pip install -e .`.

Now you can use COMET on the command line (run `comet -h` for help) and you can import functions from the COMET library into your own project. 

## Usage

### Input file

COMET needs an input file in `tsv` or `csv` format. An example input file can be seen in `data/example.csv`. 

A header is optional and the presence of a header in the input file can be flagged with `--header` when using COMET on the command line.

There is no limit on the number of molecule groups in the input file. Individual names are also optional but the column should still be present as an empty column.

### Running COMET

COMET requires setting the path to the input file and the output file, with:
* `-i`/`--input`: input file path. 
* `-o`/`--output`: output file path. 

There are several optional parameters:
* `-a`/`--alpha`: significance threshold (default: `0.05`).
* `-n`/`--num_bits`: number of bits in barcodes to construct (default: `2048`).
* `-r`/`--radius`: radius around atoms for constructing barcodes (default: `3`).
* `-c`/`--mtc`: multiple testing correction type (default: `bonferroni`).

COMET can be run on the example input as follows to get the same results as in `data/example.out`:

`comet -i ./data/example.csv --header -o ./data/example.out -a 0.1`
