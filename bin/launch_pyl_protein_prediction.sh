#!/usr/bin/env bash

echo "Activate conda environment"
source activate pyl_protein_prediction
echo ""

echo "Launch"
pythow src/gui/menu.py
echo ""
