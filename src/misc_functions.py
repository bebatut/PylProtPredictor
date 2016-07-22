#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-
# -*- coding: utf-8 -*-
from Bio import SeqIO 

def isfloat(value):
    try:
        float(value)
        return True
    except ValueError:
        return False

def isfasta(file):
	fasta = False
	for record in SeqIO.parse(file,"fasta"):
		fasta = True
	return fasta
