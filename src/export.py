#!/usr/bin/env python

import pandas as pd
from Bio import SeqIO


def export_csv(dictionary, csv_filepath, col):
    """Export a dictionary into a csv file with specific columns

    :param dictionary: dictionary with sequence ids as keys and start, end, strand, and origin sequences as values
    :param csv_filepath: path to CSV file to export
    :param col: list of columns (keys of dictionary's subdict) to conserve
    """
    dict_df = pd.DataFrame(data=dictionary).transpose()
    dict_df = dict_df[col]
    dict_df.to_csv(csv_filepath)


def export_fasta(sequences, output_filepath):
    """Export a list of SeqRecord into a fasta file
    
    :param sequences: list of SeqRecord
    :param output_filepath: path to a fasta file
    """
    with open(output_filepath, "w") as output_file:
        SeqIO.write(sequences, output_file, "fasta")