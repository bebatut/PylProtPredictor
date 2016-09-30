#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-
# -*- coding: utf-8 -*-
import os


def launch_prodigal(genome_filepath, predicted_cds_filepath):
    dirpath = os.path.dirname(predicted_cds_filepath)
    info_filepath = dirpath + "/cds_prediction_info.txt"
    log_filepath = dirpath + "/cds_prediction_log.txt"

    cmd = "prodigal"
    cmd += " -i " + genome_filepath
    cmd += " -d " + predicted_cds_filepath
    cmd += " -f gbk"
    cmd += " -g 11"
    cmd += " -o " + info_filepath
    cmd += " 2> " + log_filepath
    os.system(cmd)

def predict_cds(genome_filepath, predicted_cds_filepath,
cds_prediction_tool = "prodigal"):
    if cds_prediction_tool == "prodigal":
        launch_prodigal(genome_filepath, predicted_cds_filepath)
