#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-
# -*- coding: utf-8 -*-

import os
import sys
import re
from Bio.Blast import NCBIXML
from Bio import SeqIO

def format_ref_db(ref_db_filepath, log_file):
    ref_db_path = '.'.join(ref_db_filepath.split(".")[:-1])
    if os.path.exists(ref_db_path + '.nal'):
        log_file.write("Blast db for reference database already built!")
    else:
        log_file.write("Build Blast db for reference database\n")
        log_file.write("-------------------------------------\n")
        cmd = "makeblastdb "
        cmd += "-in " + ref_db_filepath
        cmd += " -dbtype 'prot'"
        #os.system(cmd)

def check_potential_pyl_protein(pot_pyl_prot_filepath, ref_db_filepath,
    output_dirpath, evalue):
    pot_pyl_prot_filename = os.path.basename(pot_pyl_prot_filepath)
    pot_pyl_prot_name = pot_pyl_prot_filename.split(".")[0]

    log_filepath = output_dirpath + "/" + pot_pyl_prot_name
    log_filepath += "_Pyl_protein_checking_log.txt"
    log_file = open(log_filepath, 'w')

    format_ref_db(ref_db_filepath, log_file)
