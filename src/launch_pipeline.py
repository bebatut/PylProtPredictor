#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-
# -*- coding: utf-8 -*-

from gooey import Gooey, GooeyParser
import message
import os

@Gooey(program_name="PYL protein prediction",required_cols=1)
def main():
    desc = "Predict PYL protein from an assembled or scaffold genome"
    file_help_msg = "Name of the file you want to process"

    my_cool_parser = GooeyParser(description=desc)

    my_cool_parser.add_argument("genome",
        metavar="Genome",
        help="FASTA file with an assembled or scaffold genome",
        widget="FileChooser")
    my_cool_parser.add_argument("output",
        metavar="Output directory",
        help="For all generated files (results and logs)",
        widget="DirChooser")

    my_cool_parser.add_argument('--to_launch',
        metavar="Part of whole pipelin to launch",
        default="Whole pipeline",
        choices=['Whole pipeline', 'CDS prediction',
        'PYL protein prediction', 'PYL protein checking'],
        help='')

    my_cool_parser.add_argument("--ref_db",
        metavar="Reference database",
        help="FASTA file with an reference database used to confirm potential PYL proteins", default= os.getcwd() + "data/uniref90.fasta",
        widget="FileChooser")

    #my_cool_parser.add_argument("--predicted_cds",
    #    metavar="Predicted CDS",
    #    help="FASTA file with CDS predicted using the current tool",
    #    widget="FileChooser")

    #my_cool_parser.add_argument("--pot_pyl_prot",
    #    metavar="Sequences of a predicted PYL protein",
    #    help="FASTA file with sequences of a PYL protein predicted using the current tool",
    #    widget="FileChooser")
    #my_cool_parser.add_argument('--similarity_search_tool',
    #    metavar="Similarity search tool to use",
    #    default="Diamond",
    #    choices=['Diamond', 'BLAST'],
    #    help='')

    args = my_cool_parser.parse_args()


if __name__ == '__main__':
    main()
