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

```
$ git clone git@gitlab.com:bebatut/pyl_protein_prediction.git
```

- Move into the directory

```
$ cd pyl_protein_prediction
```

- Prepare the environment

```
$ conda env create --name PylProtPredictor --file environment.yml
$ source activate PylProtPredictor
```

> To exit the environment, you can execute
> ```
> $ source deactivate
> ```
> But don't do that before running the analysis.

# Usage

```
$ PylProtPredictor --genome FILE --output PATH [options]
```

The first run will be long: the reference database should be downloaded and prepare for the similarity search.

## Database setup

If you already have the Uniref90 database on your machine, you can simply symlink it.

``` bash
ln -s /path/to/uniref90.dmnd data/uniref90.dmnd
snakemake --cleanup-metadata data.uniref90.dmnd
```

Otherwise, the pipeline will download and format it. Make sure you have at least 25GB available for the reference database. It can take several hours, depending of your connection.

# Contributors

- Cécile Hilpert
- Bérénice Batut
- Ylana Sauvaget
- Kévin Gravouil

# Support & Bug Reports

You can file a [GitHub issue](https://github.com/bebatut/PylProtPredictor/issues).
