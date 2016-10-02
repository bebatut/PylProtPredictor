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
        cmd += " -e 0.01"
        cmd += " -f 6"
        cmd += " -b 0.5"
        cmd += " --quiet"

        msg = "Running Diamond on possible sequences for the protein against the reference database..."
        print "\t", msg
        msg += "\tcommand: " + cmd
        log_file.write(msg + "\n")

        os.system(cmd)

def check_similarity_results(sim_search_output_filepath, pot_pyl_prot_filepath,
log_file, output_dirpath):
    pot_pyl_prot = {}
    for record in SeqIO.parse(pot_pyl_prot_filepath, "fasta"):
        pot_pyl_prot[record.id] = record

    conserved_seq_info = {"evalue":10, "id":0}
    with open(sim_search_output_filepath, "r") as sim_search_output_file:
        for line in sim_search_output_file.readlines():
            split_line = line[:-1].split("\t")
            evalue = float(split_line[10])
            if evalue < conserved_seq_info["evalue"]:
                conserved_seq_info["evalue"] = evalue
                conserved_seq_info["id"] = split_line[0]

    rejected_seq = []
    conserved_seq = []
    for seq_id in pot_pyl_prot:
        if seq_id == conserved_seq_info["id"]:
            conserved_seq.append(pot_pyl_prot[seq_id])
        else:
            rejected_seq.append(pot_pyl_prot[seq_id])
    rejected_seq_file = open(output_dirpath + "/rejected_protein_sequences.fasta" , "w")
    SeqIO.write(rejected_seq, rejected_seq_file, "fasta")
    rejected_seq_file.close()

    conserved_seq_file = open(output_dirpath + "/protein_sequence.fasta" , "w")
    SeqIO.write(conserved_seq, conserved_seq_file, "fasta")
    conserved_seq_file.close()

    if conserved_seq_info["id"] == 0:
        msg = "-- None of the possible sequences match the reference database --\n"
        print "\t", msg
        log_file.write(msg + "\n")
        msg = "    We can not say if the protein is able to use PYL amino acid\n"
        print "\t", msg
        log_file.write(msg + "\n")
    elif conserved_seq_info["id"].endswith("_1"):
        msg = "-- This protein does not use PYL amino acid --\n"
        print "\t", msg
        log_file.write(msg + "\n")
    else:
        msg = ">>>> This protein contain PYL amino acid <<<<"
        print "\t", msg
        log_file.write(msg + "\n")

        new_prot_description = pot_pyl_prot[conserved_seq_info["id"]].description
        msg = "\talternative protein: " + new_prot_description + "\n"
        print "\t", msg
        log_file.write(msg + "\n")

def check_potential_pyl_protein(pot_pyl_prot_filepath, ref_db_filepath,
    output_dirpath, similarity_search_tool = "diamond"):
    log_filepath = output_dirpath + "/pyl_checking_log.txt"

    with open(log_filepath, 'w') as log_file:
        format_ref_db(ref_db_filepath, log_file, similarity_search_tool)

        sim_search_output_filepath = output_dirpath + "/similarity_search.txt"
        run_similarity_search(pot_pyl_prot_filepath, ref_db_filepath,
        sim_search_output_filepath, log_file, similarity_search_tool)

        check_similarity_results(sim_search_output_filepath,
        pot_pyl_prot_filepath, log_file, output_dirpath)
