#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-
# -*- coding: utf-8 -*-

from gooey import Gooey, GooeyParser
import message
import os

@Gooey(program_name="PYL protein prediction",required_cols=1)
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
    whole_pipeline_parser.add_argument("--ref_db",
        metavar="Reference database",
        help="FASTA file with an reference database used to confirm potential PYL proteins", default= os.getcwd() + "data/uniref90.fasta",
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
    pyl_checking_parser.add_argument("output",
        metavar="Output directory",
        help="For all generated files (results and logs)",
        widget="DirChooser")
    pyl_checking_parser.add_argument('--similarity_search_tool',
        metavar="Similarity search tool to use",
        default="Diamond",
        choices=['Diamond', 'BLAST'],
        help='')

    args = parser.parse_args()
    launch_pipeline(args)


if __name__ == '__main__':
    main()
