Detection of pyrrolysine proteins
=================================

# Introduction

Pyrrolysine is an amino acid that is used in the biosynthesis of proteins in some methanogenic archaea and bacterium. It is encoded in mRNA by the UAG codon, which in most organisms is the 'amber' stop codon.

Some methanogenic archaea and bacterium have the pylT gene, which encodes an unusual transfer RNA (tRNA) with a CUA anticodon, and the pylS gene, which encodes a class II aminoacyl-tRNA synthetase that charges the pylT-derived tRNA with pyrrolysine. In some proteins, the UAG codon can then code for pyrrolysine, and no more for a stop codon.

These proteins are difficult to identify. Indeed, in CDS prediction, UAG codons are seen as stop codons. The predicted CDS are then cut when the first UAG codon is found.

Here, we propose a solution to detect pyrrolysine proteins. It relies on:

1. CDS prediction (using Prodigal)
2. Detection of possible pyrrolysine proteins in predicted CDS
    1. Identification of predicted CDS ending with a TAG codon
    2. Elongation of the CDS toward next STOP codon
    3. Conservation of all possible sequences (without or with PYL amino acid)
3. Checking each potential PYL protein
    1. Similarity search of possible sequences (without or with PYL amino acid) against known proteins
    2. Checking similarity search result: if a sequence using PYL (instead of STOP codon) has a smaller e-value, the protein is considered as a protein using PYL

# Installation

- Install the requirements
  - `git`
  - `conda`
- Clone this repository (or get the release)

```
$ git clone git@gitlab.com:bebatut/pyl_protein_prediction.git
```

- Move into the directory

```
$ cd pyl_protein_prediction
```

- Prepare environment (create `conda` environment, download database)

```
$ ./bin/prepare_environment
```

# Usage

- Launch program

```
$ ./bin/launch_pyl_protein_prediction
```

# Contributers

- Cécile Hilpert
- Bérénice Batut

# Support & Bug Reports

You can file an GitHub issue.
