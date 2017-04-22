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

Install the requirements
- `git`
- `conda`

Make sure you have at least 20GB available for the reference database (UniRef90)

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

- Download and prepare the reference database

```
$ snakemake data/uniref90.dmnd
``` 

# Usage

- Change the path to the genome in the `config.yml`
- Launch the workflow

```
$ snakemake
```

# Contributors

- Cécile Hilpert
- Bérénice Batut
- Ylana Sauvaget
- Kévin Gravouil

# Support & Bug Reports

You can file an GitHub issue.
