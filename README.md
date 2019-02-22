Detection of pyrrolysine proteins
=================================

[![CircleCI](https://circleci.com/gh/bebatut/PylProtPredictor/tree/master.svg?style=svg)](https://circleci.com/gh/bebatut/PylProtPredictor/tree/master)
[![codecov](https://codecov.io/gh/bebatut/PylProtPredictor/branch/master/graph/badge.svg?token=6KyTn6n8Bp)](https://codecov.io/gh/bebatut/PylProtPredictor)
[![DOI](https://zenodo.org/badge/88719042.svg)](https://zenodo.org/badge/latestdoi/88719042)



# Context

Pyrrolysine is an amino acid that is used in the biosynthesis of proteins in some methanogenic archaea and bacterium. It is encoded in mRNA by the UAG codon, which in most organisms is the 'amber' stop codon.

Some methanogenic archaea and bacterium have the pylT gene, which encodes an unusual transfer RNA (tRNA) with a CUA anticodon, and the pylS gene, which encodes a class II aminoacyl-tRNA synthetase that charges the pylT-derived tRNA with pyrrolysine. In some proteins, the UAG codon can then code for pyrrolysine, and no more for a STOP codon.

These proteins are difficult to identify. Indeed, in CDS prediction, UAG codons are seen as STOP codons. The predicted CDS are then cut when the first UAG codon is found.

Here, we propose a solution to detect proteins using Pyrrolisine amino acid.
Have a look to [the scheme explaining how the tool is working](doc/img/main_scheme.png).


# Installation

## Requirements

The following software are required:
- [`git`](https://git-scm.com/book/fr/v1/D%C3%A9marrage-rapide-Installation-de-Git#Installation-sur-Linux)
- [`conda`](https://conda.io/miniconda.html):

    ```
    $ make install-conda
    $ make configure-conda
    ```

## Install the tool

- Clone this repository (or get the release)

    ```
    $ git clone https://github.com/bebatut/PylProtPredictor.git
    ```

- Move into the directory

    ```
    $ cd pyl_protein_prediction
    ```

- Prepare the environment (only once)

    ```
    $ make create-env
    ```

- Activate the conda environment

    ```
    $ source activate PylProtPredictor
    ```

- Build the package

    ```
    $ make init
    ```

# Usage

```
$ source activate PylProtPredictor # once to activate the conda environment
$ pylprotpredictor --help
usage: pylprotpredictor [-h] --genome GENOME --output OUTPUT
                        [--reference_fasta_db REFERENCE_FASTA_DB]
                        [--reference_dmnd_db REFERENCE_DMND_DB]

PylProtPredictor Pipeline

optional arguments:
  -h, --help            show this help message and exit
  --genome GENOME       path to a FASTA file with full or contig sequences of
                        a genome to analyze
  --output OUTPUT       path to the output directory
  --reference_fasta_db REFERENCE_FASTA_DB
                        path to FASTA file with reference database
  --reference_dmnd_db REFERENCE_DMND_DB
                        path to Diamond formatted file with reference database
```

> To exit the environment, you can execute
> ```
> $ source deactivate
> ```
> But don't do that before running the analysis.

## Database setup

The first run will be long: the reference database should be downloaded and prepare for the similarity search.

If you already have the Uniref90 database on your machine, you can simply link it when running the main script.

Otherwise, the pipeline will download and format it. Make sure you have at least 25GB available for the reference database. It can take several hours, depending on your connection.

# Support & Bug Reports

You can file a [GitHub issue](https://github.com/bebatut/PylProtPredictor/issues).

# Contributing

First off, thanks for taking the time to contribute!

## Tests

The code is covered by tests. They are run automatically on CircleCI but we also recommend to run them locally before pushing to GitHub with:

```
$ make test
```

Any added code should be covered by new tests.

## Documentation

Documentation about ENASearch is available online at http://bebatut.fr/PylProtPredictor

To update it:

- Make the changes in `src/docs`
- Generate the doc:

    ```
    $ make doc
    ```

- Check it by opening the `docs/index.html` file in a web browser
- Propose the changes via a Pull Request

## Contributors

- Bérénice Batut
- Jean-François Brugère
- Kévin Gravouil
- Cécile Hilpert
- Ylana Sauvaget

# Citation

You can cite the latest release on [Zenodo](https://zenodo.org/record/2575708) 

