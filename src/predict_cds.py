#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-
# -*- coding: utf-8 -*-
import os


def launch_prodigal(genome_filepath, result_dirpath):
    genome_filename = os.path.basename(genome_filepath)
    genome_name = genome_filename.split(".")[0]

    result_filepath = result_dirpath + "/" + genome_name
    result_filepath += "_predicted_CDS.fasta"

    info_filepath = result_dirpath + "/" + genome_name
    info_filepath += "_CDS_prediction_info.txt"

    log_filepath = result_dirpath + "/" + genome_name
    log_filepath += "_CDS_prediction_log.txt"

    cmd = "prodigal"
    cmd += " -i " + genome_filepath
    cmd += " -d " + result_filepath
    cmd += " > " + info_filepath
    cmd += " 2> " + log_filepath
    os.system(cmd)
