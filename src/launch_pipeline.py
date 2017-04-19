#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-
# -*- coding: utf-8 -*-

from gooey import Gooey, GooeyParser
import message
import os
import misc_functions
import predict_cds
import predict_pyl_proteins
import check_pyl_protein

def launch_cds_prediction(genome_filepath, output_dir):
    print "Predict CDS..."
    if not misc_functions.isfasta(genome_filepath):
        raise ValueError("The file with the genome is not a FASTA file")

    predicted_cds_filepath = output_dir + "/predicted_cds.fasta"
    predict_cds.predict_cds(genome_filepath, predicted_cds_filepath,
    "prodigal")
    return predicted_cds_filepath

def launch_pyl_protein_prediction(genome_filepath, predicted_cds_filepath,
output_dir):
    print "Predict PYL protein..."
    if not misc_functions.isfasta(genome_filepath):
        raise ValueError("The file with the genome is not a FASTA file")
    if not misc_functions.isfasta(predicted_cds_filepath):
        raise ValueError("The file with the predicted CDS is not a FASTA file")

    predict_pyl_proteins.predict_pyl_proteins(genome_filepath,
    predicted_cds_filepath, output_dir)

def launch_pyl_protein_checking(pot_pyl_prot_filepath, ref_db_filepath,
similarity_search_tool, output_dir):
    if not misc_functions.isfasta(pot_pyl_prot_filepath):
        raise ValueError("The file with the genome is not a FASTA file")

    prot_name = os.path.basename(output_dir)
    print "Checking if " + prot_name + " is a PYL protein (it can take more than 1 hour)..."

    check_pyl_protein.check_potential_pyl_protein(pot_pyl_prot_filepath, ref_db_filepath, output_dir, similarity_search_tool)

def run_whole_pipeline(genome_filepath, output_dir, ref_db_filepath):
    predicted_cds_filepath = launch_cds_prediction(genome_filepath, output_dir)

    pyl_protein_dir = output_dir + "/potential_pyl_prot"
    if not os.path.exists(pyl_protein_dir):
        os.mkdir(pyl_protein_dir)
    launch_pyl_protein_prediction(genome_filepath, predicted_cds_filepath,
    pyl_protein_dir)

    listing = next(os.walk(pyl_protein_dir))[1]
    for subdir in listing:
        subdirpath = pyl_protein_dir + "/" + subdir
        pot_pyl_prot_filepath = subdirpath + "/potential_sequences.fasta"
        launch_pyl_protein_checking(pot_pyl_prot_filepath, ref_db_filepath,
            "diamond", subdirpath)

def launch_pipeline(args):
    if args.subparser_name == "run-whole-pipeline":
        run_whole_pipeline(args.genome, args.output, args.ref_db)
    elif args.subparser_name == "predict-cds":
        predict_cds(args.genome, args.output)
    elif args.subparser_name == "predict-pyl-proteins":
        predict_pyl_proteins(args.genome, args.predicted_cds, args.output)
    elif args.subparser_name == "check-pyl-protein":
        check_pyl_protein(args.pot_pyl_prot, args.ref_db,
        args.similarity_search_tool, args.output)
    else:
        raise ValueError("Wrong command")

@Gooey(
    optional_cols=2,
    program_name="PYL protein prediction",
    default_size=(900, 600)
)
def main():
    desc = "Predict PYL proteins from an assembled or scaffold genome"
    file_help_msg = "Name of the file you want to process"

    parser = GooeyParser(description=desc)

    subs = parser.add_subparsers(dest='subparser_name')

    whole_pipeline_parser = subs.add_parser("run-whole-pipeline",
        help='Whole pipeline will be executed with CDS prediction, PYL protein prediction and checking')
    whole_pipeline_parser.add_argument("genome",
        metavar="Genome",
        help="FASTA file with an assembled or scaffold genome",
        widget="FileChooser")
    whole_pipeline_parser.add_argument("output",
        metavar="Output directory",
        help="For all generated files (results and logs)",
        widget="DirChooser")
    whole_pipeline_parser.add_argument("ref_db",
        metavar="Reference database",
        help="FASTA file with an reference database used to confirm potential PYL proteins", default= os.getcwd() + "/data/uniref90.fasta",
        widget="FileChooser")

    cds_prediction_parser = subs.add_parser("predict-cds",
        help='CDS will predicted for the genome using Prodigal')
    cds_prediction_parser.add_argument("genome",
        metavar="Genome",
        help="FASTA file with an assembled or scaffold genome",
        widget="FileChooser")
    cds_prediction_parser.add_argument("output",
        metavar="Output directory",
        help="For all generated files (results and logs)",
        widget="DirChooser")

    pyl_prediction_parser = subs.add_parser("predict-pyl-proteins",
        help='PYL proteins will be predicted based on predicted CDS. The TAG ending proteins will be extended to next STOP codons. One FASTA file will be generated for each possible protein with the sequences to check to determine if the protein use PYL instead of STOP for TAG codons.')
    pyl_prediction_parser.add_argument("predicted_cds",
        metavar="Predicted CDS",
        help="FASTA file with CDS predicted using the current tool",
        widget="FileChooser")
    pyl_prediction_parser.add_argument("genome",
        metavar="Genome",
        help="FASTA file with an assembled or scaffold genome",
        widget="FileChooser")
    pyl_prediction_parser.add_argument("output",
        metavar="Output directory",
        help="For all generated files (results and logs)",
        widget="DirChooser")

    pyl_checking_parser = subs.add_parser("check-pyl-protein",
        help='A potential PYL protein will be checked. A similarity search is runned for the possible sequences of the protein against a reference database. If an extended sequence (i.e. where a TAG codon was transformed in PYL) has better similarity search results against the reference database, the protein is considered as a PYL protein. Otherwise, the protein is not protein using the PYL amino acid.')
    pyl_checking_parser.add_argument("pot_pyl_prot",
        metavar="Sequences of a predicted PYL protein",
        help="FASTA file with sequences of a PYL protein predicted using the current tool",
        widget="FileChooser")
    pyl_checking_parser.add_argument('similarity_search_tool',
        metavar="Similarity search tool to use",
        default="Diamond",
        choices=['Diamond', 'BLAST'],
        help='')
    pyl_checking_parser.add_argument("ref_db",
        metavar="Reference database",
        help="FASTA file with an reference database used to confirm potential PYL proteins", default= os.getcwd() + "/data/uniref90.fasta",
        widget="FileChooser")
    pyl_checking_parser.add_argument("output",
        metavar="Output directory",
        help="For all generated files (results and logs)",
        widget="DirChooser")


    args = parser.parse_args()
    launch_pipeline(args)


if __name__ == '__main__':
    main()
