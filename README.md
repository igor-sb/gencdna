# gencDNA library toolkit

[![CI](https://github.com/igor-sb/gencdna/actions/workflows/ci.yml/badge.svg)](https://github.com/igor-sb/gencdna/actions)
[![codecov](https://codecov.io/gh/igor-sb/gencdna/branch/PD-17-initial-setup/graph/badge.svg?token=AJ97Y57FVK)](https://codecov.io/gh/igor-sb/gencdna)

Code for detecting genomic cDNAs (gencDNAs) from various WGS or PCR-based DNA
sequencing data. Many of these scripts are ran in a pipeline together with
other bioinformatic tools such as [samtools](https://www.htslib.org/), 
[bedtools](https://bedtools.readthedocs.io/en/latest/index.html) and 
[bwa](https://github.com/lh3/bwa). So, these
may need to be installed if you want to run the entire pipeline.

## Documentation

See [Documentation](https://igor-sb.github.io/gencdna/) for detailed instructions
how to run these scripts.

## Installation

Install [Poetry](https://python-poetry.org/) then run:

```
git clone https://github.com/igor-sb/gencdna.git
cd gencdna
make install
```

This will create a new virtual environment that you can activate or deactivate
using `poetry shell` when inside the `gencdna` folder.