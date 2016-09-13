#!/usr/bin/env python
# -*- coding: utf-8 -*-

# from Bio import SeqIO
# from Bio import SeqRecord
# from Bio.Seq import Seq
# from Bio.SeqRecord import SeqRecord
# from Bio.Alphabet import generic_protein
import re
import os
import sys
from Bio import SeqIO
from Bio.Data import CodonTable
import misc_functions

def extract_predicted_cds_position(predicted_cds_filepath):
    pred_cds_pos = {}
    pred_cds_pos['forward'] = {}
    pred_cds_pos['forward']['details'] = {}
    pred_cds_pos['forward']['order'] = []
    pred_cds_pos['reverse'] = {}
    pred_cds_pos['reverse']['details'] = {}
    pred_cds_pos['reverse']['order'] = []
    pred_cds_nb = 0
    for record in SeqIO.parse(predicted_cds_filepath,"fasta"):
        pred_cds_nb += 1
        split_description = record.description.split("#")
        seq_id = split_description[0][:-1]
        start = split_description[1].replace(" ","")
        end = split_description[2].replace(" ","")
        strand = split_description[3].replace(" ","")
        if strand == "-1":
            pred_cds_pos['reverse']['details'].setdefault(seq_id, {})
            pred_cds_pos['reverse']['details'][seq_id]['start'] = int(start)
            pred_cds_pos['reverse']['details'][seq_id]['end'] = int(end)
            pred_cds_pos['reverse']['details'][seq_id]['seq'] = record.seq
            pred_cds_pos['reverse']['order'].append(seq_id)
        else:
            pred_cds_pos['forward']['details'].setdefault(seq_id, {})
            pred_cds_pos['forward']['details'][seq_id]['start'] = int(start)
            pred_cds_pos['forward']['details'][seq_id]['end'] = int(end)
            pred_cds_pos['forward']['details'][seq_id]['seq'] = record.seq
            pred_cds_pos['forward']['order'].append(seq_id)
    return pred_cds_pos, pred_cds_nb

def identify_tag_ending_protein(pred_cds_pos):
    tag_ending_protein = {}
    tag_ending_protein['forward'] = {}
    tag_ending_protein['reverse'] = {}
    tag_ending_prot_nb = 0
    for strand in pred_cds_pos:
        for cds in pred_cds_pos[strand]['details']:
            seq = pred_cds_pos[strand]['details'][cds]['seq']
            if str(seq).endswith("TAG"):
                tag_ending_prot_nb += 1
                tag_ending_protein[strand].setdefault(cds, {})
                start = pred_cds_pos[strand]['details'][cds]['start']
                tag_ending_protein[strand][cds]["start"] = start
                end = pred_cds_pos[strand]['details'][cds]['end']
                tag_ending_protein[strand][cds]["ends"] = [end]
                seq = pred_cds_pos[strand]['details'][cds]['seq']
                tag_ending_protein[strand][cds]["seqs"] = [seq]
    return tag_ending_protein, tag_ending_prot_nb

def extend_to_next_stop_codon(current_end, genome):
    genome = str(genome.seq)
    genome_size = len(genome)
    stop_codons = CodonTable.unambiguous_dna_by_id[1].stop_codons

    new_end = current_end
    to_continue = (new_end+3) < genome_size
    new_ends = []
    #print 'last codon:'
    #print genome[(new_end-3):new_end]
    while to_continue:
        codon = str(genome[new_end:(new_end + 3)])
        #print "\t",codon
        if codon not in stop_codons :
            new_end += 3
            to_continue = (new_end+3) < genome_size
        else:
            new_end += 3
            new_ends.append(new_end)
            if codon != 'TAG':
                to_continue = False
    return new_ends

def extend_tag_ending_proteins(tag_ending_protein, pred_cds_pos,
    genome_filepath):
    genome = SeqIO.read(genome_filepath, "fasta")
    genome_size = len(genome.seq)

    strand = "forward"
    cds_nb = len(pred_cds_pos[strand]['order'])
    for cds in tag_ending_protein[strand]:
        start = tag_ending_protein[strand][cds]['start']
        end = tag_ending_protein[strand][cds]['ends'][0]
        new_ends = extend_to_next_stop_codon(end, genome)
        tag_ending_protein[strand][cds]['ends'] += new_ends
        for new_end in new_ends:
            new_seq = genome.seq[start:new_end]
            tag_ending_protein[strand][cds]['seqs'].append(new_seq)

    strand = "reverse"
    reverse_complement_genome = genome.reverse_complement()
    cds_nb = len(pred_cds_pos[strand]['order'])
    for cds in tag_ending_protein[strand]:
        start = tag_ending_protein[strand][cds]['start']
        end = tag_ending_protein[strand][cds]['ends'][0]
        reverse_start = (genome_size-end)
        reverse_end = (genome_size-start+1)
        new_ends = extend_to_next_stop_codon(reverse_end,
            reverse_complement_genome)
        for new_end in new_ends:
            new_seq = reverse_complement_genome.seq[reverse_start:new_end]
            tag_ending_protein[strand][cds]['seqs'].append(new_seq)

    return tag_ending_protein

def predict_pyl_protein_on_scaffold_genome(genome_filepath,
    predicted_cds_filepath, output_dirpath, log_file):
    print "launch_pyl_protein_from_scaffold"


def predict_pyl_protein_on_assembled_genome(genome_filepath,
    predicted_cds_filepath, output_dirpath, log_file):
    pred_cds_pos, pred_cds_nb = extract_predicted_cds_position(
    predicted_cds_filepath)
    log_file.write("Number of predicted CDS: ")
    log_file.write(str(pred_cds_nb) + "\n")

    tag_ending_protein, tag_ending_prot_nb = identify_tag_ending_protein(
    pred_cds_pos)
    log_file.write("Number of TAG-ending predicted CDS: ")
    log_file.write(str(tag_ending_prot_nb) + "\n")

    tag_ending_protein = extend_tag_ending_proteins(tag_ending_protein,
    pred_cds_pos, genome_filepath)

def predict_pyl_proteins(genome_filepath, predicted_cds_filepath,
    output_dirpath):
    genome_filename = os.path.basename(genome_filepath)
    genome_name = genome_filename.split(".")[0]

    log_filepath = output_dirpath + "/" + genome_name
    log_filepath += "_Pyl_protein_prediction_log.txt"
    log_file = open(log_filepath, 'w')

    if misc_functions.isscaffold(genome_filepath):
        log_file.write("Pyl prediction on scaffold genome\n")
        log_file.write("---------------------------------\n")
        predict_pyl_protein_on_scaffold_genome(genome_filepath,
        predicted_cds_filepath, output_dirpath, log_file)
    else:
        log_file.write("Pyl prediction on assembled genome\n")
        log_file.write("----------------------------------\n")
        predict_pyl_protein_on_assembled_genome(genome_filepath,
        predicted_cds_filepath, output_dirpath, log_file)

    log_file.close()

predict_pyl_protein_with_assemble_genome("../data/MalvusMx1201_GENOME_FASTA.fsa",
"../data/MalvusMx1201_GENOME_FASTA_predicted_CDS.fasta", "../data")
