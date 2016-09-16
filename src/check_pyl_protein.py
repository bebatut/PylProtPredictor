#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-
# -*- coding: utf-8 -*-

import os
import sys
from Bio import SeqIO

def format_ref_db(ref_db_filepath, log_file, similarity_search_tool):
    if similarity_search_tool == "blast":
        if os.path.exists(ref_db_filepath + '.pal'):
            msg = "Blast db for reference database is already built!"
            log_file.write(msg + "\n")
            print "\t", msg
        else:
            cmd = "makeblastdb "
            cmd += "-in " + ref_db_filepath
            cmd += " -dbtype 'prot'"

            msg = "Building Blast db for reference database...\n"
            print "\t", msg
            msg += "\tcommand: " + cmd
            log_file.write(msg + "\n")

            os.system(cmd)
    else:
        ref_db_path = ".".join(ref_db_filepath.split(".")[:-1])
        if os.path.exists(ref_db_path + '.dmnd'):
            msg = "Diamond db for reference database is already built!"
            log_file.write(msg + "\n")
            print "\t", msg
        else:
            cmd = "diamond makedb"
            cmd += " --in " + ref_db_filepath
            cmd += " --db " + ref_db_path
            cmd += " -b 0.5"
            cmd += " --quiet"

            msg = "Building Diamond db for reference database...\n"
            print "\t", msg
            msg += "\tcommand: " + cmd
            log_file.write(msg + "\n")

            os.system(cmd)

def run_similarity_search(pot_pyl_prot_filepath, ref_db_filepath,
    output_filepath, log_file, similarity_search_tool):
    if similarity_search_tool == "blast":
        cmd = "blastp "
        cmd += " -db " + ref_db_filepath
        cmd += " -query " + pot_pyl_prot_filepath
        cmd += " -out " + output_filepath
        cmd += " -evalue 1e-10"
        cmd += " -outfmt 6"
        cmd += " -max_target_seqs 1"
        cmd += " -task blastp-fast"

        msg = "Running Blast on possible sequences for the protein against the reference database..."
        print "\t", msg
        msg += "\tcommand: " + cmd
        log_file.write(msg + "\n")

        os.system(cmd)
    else:
        ref_db_path = ".".join(ref_db_filepath.split(".")[:-1])
        cmd = "diamond blastp"
        cmd += " -d " + ref_db_path
        cmd += " -q " + pot_pyl_prot_filepath
        cmd += " -o " + output_filepath
        cmd += " -k 1"
        cmd += " -e 1e-10"
        cmd += " -f 6"
        cmd += " -b 0.5"
        cmd += " --quiet"

        msg = "Running Diamond on possible sequences for the protein against the reference database..."
        print "\t", msg
        msg += "\tcommand: " + cmd
        log_file.write(msg + "\n")

        os.system(cmd)

def check_similarity_results(sim_search_output_filepath, pot_pyl_prot_filepath,
log_file):
    with open(sim_search_output_filepath, "r") as sim_search_output_file:
        smallest_evalue = {"value":10, "id":0}
        for line in sim_search_output_file.readlines():
            split_line = line[:-1].split("\t")
            evalue = float(split_line[10])
            if evalue < smallest_evalue["value"]:
                smallest_evalue["value"] = evalue
                smallest_evalue["id"] = split_line[0]

    if smallest_evalue["id"].endswith("_1"):
        msg = "-- This protein does not use PYL as amino acid --\n"
        print "\t", msg
        log_file.write(msg + "\n")
    else:
        msg = ">>>> This protein seems to use PYL as amino acid <<<<"
        print "\t", msg
        log_file.write(msg + "\n")

        for record in SeqIO.parse(pot_pyl_prot_filepath, "fasta"):
            if record.id == smallest_evalue["id"]:
                new_prot_description = record.description

        msg = "\talternative protein: " + new_prot_description + "\n"
        print "\t", msg
        log_file.write(msg + "\n")

def check_potential_pyl_protein(pot_pyl_prot_filepath, ref_db_filepath,
    output_dirpath, similarity_search_tool = "diamond"):

    pot_pyl_prot_filebase = os.path.basename(pot_pyl_prot_filepath)
    prot_name = pot_pyl_prot_filebase.split(".")[0]

    log_filepath = output_dirpath + "/" +  prot_name + "_pyl_checking_log.txt"

    with open(log_filepath, 'w') as log_file:
        format_ref_db(ref_db_filepath, log_file, similarity_search_tool)

        sim_search_output_filepath = output_dirpath + "/"
        sim_search_output_filepath += prot_name + "_similarity_search.txt"
        run_similarity_search(pot_pyl_prot_filepath, ref_db_filepath,
        sim_search_output_filepath, log_file, similarity_search_tool)

        check_similarity_results(sim_search_output_filepath,
        pot_pyl_prot_filepath, log_file)
