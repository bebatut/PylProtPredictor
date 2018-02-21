import json

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


def export_json(in_dict, output_filepath):
    """Export a dictionary into a JSON file

    :param in_dict: a dictionary
    :param output_filepath: path to the file to store the pickled object
    """
    with open(output_filepath, 'w') as f:
        f.write(json.dumps(in_dict, sort_keys=True))


def import_json(input_filepath):
    """Import a dictionary from a JSON file

    :param input_filepath: path to a JSON file

    :return: dictionary extracted from the JSON file
    """
    with open(input_filepath, 'r') as f:
        d = json.loads(f.read())
    return d
