#!/usr/bin/env bash

echo "Create conda environment"
conda env create -f pyl_protein_prediction.yml

echo "Download UniRef90 database (it can take several hours)"
wget ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz
mv uniref90.fasta.gz data/