Detection of pyrrolysine proteins
=================================


# Introduction

Pyrrolysine is an amino acid that is used in the biosynthesis of proteins in some methanogenic archaea and bacterium. It is encoded in mRNA by the UAG codon, which in most organisms is the 'amber' stop codon.

Some methanogenic archaea and bacterium have the pylT gene, which encodes an unusual transfer RNA (tRNA) with a CUA anticodon, and the pylS gene, which encodes a class II aminoacyl-tRNA synthetase that charges the pylT-derived tRNA with pyrrolysine. In some proteins, the UAG codon can then code for pyrrolysine, and no more for a STOP codon.

These proteins are difficult to identify. Indeed, in CDS prediction, UAG codons are seen as STOP codons. The predicted CDS are then cut when the first UAG codon is found.

Here, we propose a solution to detect proteins using Pyrrolisine amino acid.
Have a look to [the scheme explaining how the tool is working](doc/img/main_scheme.png).


# Installation

## Requirements

The following software are required:
- [`git`](https://git-scm.com/book/fr/v1/D%C3%A9marrage-rapide-Installation-de-Git#Installation-sur-Linux)
- [`conda`](https://conda.io/miniconda.html)

## Install the tool

- Clone this repository (or get the release)

``` bash
$ git clone https://github.com/bebatut/PylProtPredictor.git
# OR
$ wget https://github.com/bebatut/PylProtPredictor/archive/master.zip
$ unzip master.zip
```

- Move into the directory

``` bash
$ cd PylProtPredictor/
```

- Set up the conda environment

``` bash
$ conda env create --file environment.yml
```


# Run the PylProtPrediction pipeline

- Enter the conda environment

``` bash
$ source activate PylProtPredictor
```

## Snakemake introduction

[Snakemake](https://snakemake.readthedocs.io/en/stable/) is a pipeline management utility. Each pipeline step is described by a `rule`. If `ruleB` depends on `ruleA` but `ruleA` has not been executed yet, `snakemake` will automatically do it first. If `ruleA` has already been executed, `snakemake` will only execute `ruleB`. If `ruleB` has already been executed then there is nothing to do.

Here is a brief description of the pipeline.

``` bash
$ snakemake help      # print the (temporary) help
$ snakemake --list    # list all available rules (+ description)
```

## Database setup

If you already have the Uniref90 database on your machine, you can simply symlink it.

``` bash
ln -s /path/to/uniref90.dmnd data/uniref90.dmnd
snakemake --cleanup-metadata data.uniref90.dmnd
```

Otherwise, the pipeline will download and format it. Make sure you have at least 25GB available for the reference database. It can take several hours, depending of your connection. If you want to install the dabase without running the whole pipeline, run:

``` bash
$ snakemake prepare_database
```

## PylProtPredictor pipeline

- Change the path to your genome in the `config.yml`
- Launch the workflow

``` bash
$ snakemake
    --cores           # use all available CPU
    -p                # print executed commands
    PylProtPredictor  # run the complete pipeline
$ xdg-open results/YOURGENOME/report.html # open the report in your web browser
```

- Exit the conda environment
``` bash
$ source deactivate
```

## `config.yaml`

This file describes the user-defined and global parameters. You may edit it according to your needs before running the pipeline, especially the `genome`, `output_dir` and `max_threads` fields.


# Contributors

- Cécile Hilpert
- Bérénice Batut
- Ylana Sauvaget
- Kévin Gravouil

# Support & Bug Reports

You can file a [GitHub issue](https://github.com/bebatut/PylProtPredictor/issues).
