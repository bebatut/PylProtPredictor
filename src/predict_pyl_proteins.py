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
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Data import CodonTable
import misc_functions


strands = ['forward', 'reverse']

def transform_strand(strand_id):
    if strand_id == "-1":
        return strands[-1]
    else:
        return strands[0]

def extract_predicted_cds(predicted_cds_filepath):
    pred_cds = {}
    pred_cds_nb = 0
    for record in SeqIO.parse(predicted_cds_filepath,"fasta"):
        pred_cds_nb += 1
        split_description = record.description.split("#")
        seq_id = split_description[0][:-1]
        origin_seq = "_".join(seq_id.split("_")[:-1])
        start = split_description[1].replace(" ","")
        end = split_description[2].replace(" ","")
        strand = transform_strand(split_description[3].replace(" ",""))

        pred_cds.setdefault(origin_seq, {})
        pred_cds[origin_seq].setdefault(strand, {})
        pred_cds[origin_seq][strand].setdefault('details', {})
        pred_cds[origin_seq][strand].setdefault('order', [])

        pred_cds[origin_seq][strand]['details'].setdefault(seq_id, {})
        pred_cds[origin_seq][strand]['details'][seq_id]['start'] = int(start)
        pred_cds[origin_seq][strand]['details'][seq_id]['end'] = int(end)
        pred_cds[origin_seq][strand]['details'][seq_id]['seq'] = record.seq
        pred_cds[origin_seq][strand]['order'].append(seq_id)

    return pred_cds, pred_cds_nb

def identify_tag_ending_proteins(pred_cds):
    tag_ending_prot = {}
    tag_ending_prot_nb = 0
    for origin_seq in pred_cds:
        tag_ending_prot[origin_seq] = {}
        for strand in pred_cds[origin_seq]:
            tag_ending_prot[origin_seq][strand] = {}
            for i in range(len(pred_cds[origin_seq][strand]['order'])):
                cds_id = pred_cds[origin_seq][strand]['order'][i]
                seq = pred_cds[origin_seq][strand]['details'][cds_id]['seq']
                if str(seq).endswith("TAG"):
                    tag_ending_prot_nb += 1
                    start = pred_cds[origin_seq][strand]['details'][cds_id]['start']
                    end = pred_cds[origin_seq][strand]['details'][cds_id]['end']

                    tag_ending_prot[origin_seq][strand].setdefault(cds_id, {})
                    tag_ending_prot[origin_seq][strand][cds_id]["start"] = start
                    tag_ending_prot[origin_seq][strand][cds_id]["end"] = end
                    tag_ending_prot[origin_seq][strand][cds_id]["seq"] = seq
                    tag_ending_prot[origin_seq][strand][cds_id]["order_id"] = i
    return tag_ending_prot, tag_ending_prot_nb

def extract_origin_seq(genome_filepath):
    origin_seqs = {}
    for record in SeqIO.parse(genome_filepath,"fasta"):
        origin_seqs.setdefault(record.id, {})
        origin_seqs[record.id]["genome"] = record
        origin_seqs[record.id]["length"] = len(record.seq)
        origin_seqs[record.id]["rev_comp_genome"] = record.reverse_complement()
    return origin_seqs

def extend_to_next_stop_codon(current_end, genome, next_cds_end):
    genome = str(genome.seq)
    genome_size = len(genome)
    stop_codons = CodonTable.unambiguous_dna_by_id[1].stop_codons

    new_end = current_end
    to_continue = (new_end+3) < genome_size
    new_ends = []
    while to_continue:
        codon = str(genome[new_end:(new_end + 3)])
        if codon not in stop_codons :
            new_end += 3
            to_continue = (new_end+3) < genome_size
        else:
            new_end += 3
            if new_end >= next_cds_end:
                to_continue = False
            else:
                new_ends.append(new_end)
                if codon != 'TAG':
                    to_continue = False
    return new_ends

def extract_potential_pyl_proteins(tag_ending_prot, pred_cds, genome_filepath):
    origin_seqs = extract_origin_seq(genome_filepath)

    pot_pyl_prot_nb = 0
    pot_pyl_prot = {}
    for origin_seq in tag_ending_prot:
        for strand in tag_ending_prot[origin_seq]:
            pred_cds_nb = len(pred_cds[origin_seq][strand]['order'])
            for cds in tag_ending_prot[origin_seq][strand]:
                start = tag_ending_prot[origin_seq][strand][cds]['start']
                end = tag_ending_prot[origin_seq][strand][cds]['end']
                seq = tag_ending_prot[origin_seq][strand][cds]['seq']
                order_id = tag_ending_prot[origin_seq][strand][cds]['order_id']

                genome_size = origin_seqs[origin_seq]["length"]
                if strand == "reverse":
                    previous_start = start
                    start = (genome_size-end)+1
                    end = (genome_size-previous_start+1)
                    genome = origin_seqs[origin_seq]["rev_comp_genome"]
                    next_id = order_id - 1
                    if next_id >= 0:
                        next_cds_id = pred_cds[origin_seq][strand]['order'][next_id]
                        next_cds_end = pred_cds[origin_seq][strand]['details'][next_cds_id]['start']
                    else:
                        next_cds_end = genome_size
                    next_cds_end = (genome_size-next_cds_end+1)
                else:
                    genome = origin_seqs[origin_seq]["genome"]
                    next_id = order_id + 1
                    if next_id < pred_cds_nb:
                        next_cds_id = pred_cds[origin_seq][strand]['order'][next_id]
                        next_cds_end = pred_cds[origin_seq][strand]['details'][next_cds_id]['end']
                    else:
                        next_cds_end = genome_size

                new_ends = extend_to_next_stop_codon(end, genome, next_cds_end)
                if len(new_ends) > 0:
                    pot_pyl_prot_nb += 1

                    pot_pyl_prot.setdefault(cds, {})
                    pot_pyl_prot[cds]["strand"] = strand
                    pot_pyl_prot[cds]["origin_seq"] = origin_seq
                    pot_pyl_prot[cds]["potential_seq"] = []

                    if strand == "reverse":
                        pot_pyl_prot[cds]["potential_seq"].append(
                        {"start": genome_size - end + 1,
                        "end": genome_size - start + 1, "seq": seq})
                    else:
                        pot_pyl_prot[cds]["potential_seq"].append(
                        {"start": start, "end": end, "seq": seq})

                    for new_end in new_ends:
                        new_start = start
                        new_seq = genome.seq[(start-1):new_end]

                        if strand == "reverse":
                            new_start = genome_size - new_end + 1
                            new_end = genome_size - start + 1

                        pot_pyl_prot[cds]["potential_seq"].append({"start": new_start, "end": new_end, "seq": new_seq})

    return pot_pyl_prot, pot_pyl_prot_nb

def find_stop_codon_pos_in_seq(string):
    stop_codon_pos = []
    for i in range(len(string)-1):
        if string.startswith('*', i):
            stop_codon_pos.append(i)
    return stop_codon_pos

def translate(seq):
    translated_seq = seq.translate()
    str_seq = str(seq)
    n = 3
    codons = [str_seq[i:i+n] for i in range(0, len(str_seq), n)]
    for i in find_stop_codon_pos_in_seq(str(translated_seq)):
        if codons[i] != "TAG":
            raise ValueError("Stop codon found inside a sequence")
        mutable_seq = translated_seq.tomutable()
        mutable_seq[i] = "O"
        translated_seq = mutable_seq.toseq()
    return translated_seq


def save_potential_pyl_proteins(pot_pyl_prot, pyl_protein_dir):
    for prot_id in pot_pyl_prot:
        sequences = []

        count = 0

        for potential_seq in pot_pyl_prot[prot_id]["potential_seq"]:
            count += 1
            seq_id = prot_id + "_" + str(count)
            description = "# origin_seq: " + pot_pyl_prot[prot_id]["origin_seq"]
            description += " # strand: " + pot_pyl_prot[prot_id]["strand"]
            description += " # start: " + str(potential_seq["start"])
            description += " # end: " + str(potential_seq["end"])

            translated_seq = translate(potential_seq["seq"])

            sequences.append(SeqRecord(translated_seq, id=seq_id,
             description=description))

        output_file = open(pyl_protein_dir + "/" + prot_id + ".fasta" , "w")
        SeqIO.write(sequences, output_file, "fasta")
        output_file.close()
        #print

def predict_pyl_proteins(genome_filepath, predicted_cds_filepath,
    output_dirpath):
    genome_filename = os.path.basename(genome_filepath)
    genome_name = genome_filename.split(".")[0]

    pyl_protein_dir = output_dirpath + "/" + genome_name + "_potential_pyl_prot"
    if not os.path.exists(pyl_protein_dir):
        os.mkdir(pyl_protein_dir)

    log_filepath = output_dirpath + "/" + genome_name
    log_filepath += "_Pyl_protein_prediction_log.txt"
    log_file = open(log_filepath, 'w')

    pred_cds, pred_cds_nb = extract_predicted_cds(predicted_cds_filepath)
    log_file.write("Number of predicted CDS: ")
    log_file.write(str(pred_cds_nb) + "\n")

    tag_ending_prot, tag_ending_prot_nb = identify_tag_ending_proteins(
    pred_cds)
    log_file.write("Number of TAG-ending predicted CDS: ")
    log_file.write(str(tag_ending_prot_nb) + "\n")

    pot_pyl_prot, pot_pyl_prot_nb = extract_potential_pyl_proteins(
    tag_ending_prot, pred_cds, genome_filepath)
    log_file.write("Number of potential Pyl proteins:")
    log_file.write(" " + str(pot_pyl_prot_nb) + "\n")

    save_potential_pyl_proteins(pot_pyl_prot, pyl_protein_dir)

    log_file.close()

predict_pyl_proteins("../data/138T0_sorted-out_29-scaf.fna",
"../results/138T0/138T0_sorted-out_29-scaf_predicted_CDS.fasta",
"../results/138T0/")
