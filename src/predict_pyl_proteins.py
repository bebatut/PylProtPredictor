#!/usr/bin/env python

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Data import CodonTable
import pandas as pd


strands = ['forward', 'reverse']

def transform_strand(strand_id):
    """
    Transform strand from numerical value to string value
    """
    if strand_id == "-1":
        return strands[-1]
    else:
        return strands[0]


def export_csv(dictionary, csv_filepath, col):
    """
    """
    dict_df = pd.DataFrame(data=dictionary).transpose()
    dict_df = dict_df[col]
    dict_df.to_csv(csv_filepath)


def extract_predicted_cds(predicted_cds_filepath, predicted_cds_info_path):
    """
    Extract the predicted CDS
    """
    pred_cds = {}
    predicted_cds_info = {}
    for record in SeqIO.parse(predicted_cds_filepath,"fasta"):
        split_description = record.description.split("#")
        seq_id = split_description[0][:-1]
        origin_seq = "_".join(seq_id.split("_")[:-1])
        if seq_id.find("|") != -1:
            seq_id = "cds_%s" % (seq_id.split("_")[-1])
        start = split_description[1].replace(" ","")
        end = split_description[2].replace(" ","")
        strand = transform_strand(split_description[3].replace(" ",""))
        pred_cds.setdefault(origin_seq, {})
        pred_cds[origin_seq].setdefault(strand, {
            'details': {},
            'order': []})

        pred_cds[origin_seq][strand]['details'][seq_id] = {
            'start': int(start),
            'end': int(end),
            'seq': record.seq
        }
        pred_cds[origin_seq][strand]['order'].append(seq_id)
        predicted_cds_info[seq_id] = {
            "start": start,
            "end": end,
            "strand": strand,
            "origin_seq": origin_seq}

    export_csv(
        predicted_cds_info,
        predicted_cds_info_path,
        ["start", "end", "strand", "origin_seq"])
    return pred_cds, len(predicted_cds_info.keys())

def number_stop_start(predicted_cds_filepath):
    """
    Number and percentage of START and STOP codon
    """
    i = 0
    nb_stop ={"TAG":0, "TAA":0, "TGA":0} 
    nb_start={"ATG":0}
    for seq_record in SeqIO.parse("predicted_cds.fasta", "fasta"): 
        print(seq_record.id)
        i+=1
        codon_stop=seq_record.seq[-3:]
        codon_start=seq_record.seq[:3]
        print(codon_stop)
        print(codon_start)
    try: 
        nb_stop[codon_stop]+=1
        nb_start[codon_start]+=1 
    except KeyError: # Key is not present 
        pass 
print(nb_stop)
print(nb_start)
print("le fichier compte", i, "séquences.")   	
for cle, valeur in nb_stop.items():
    print("\nLe motif {} représente {}% des séquences.".format(cle, (valeur/i)*100))

for cle, valeur in nb_start.items():
    print("\nLe motif {} représente {}% des séquences.".format(cle, (valeur/i)*100))


def identify_tag_ending_proteins(pred_cds, tag_ending_cds_info_path):
    """
    Identify CDS ending with TAG STOP codon
    """
    tag_ending_prot = {}
    tag_ending_cds_info = {}
    for origin_seq in pred_cds:
        tag_ending_prot[origin_seq] = {}
        origin_cds = pred_cds[origin_seq]
        for strand in origin_cds:
            tag_ending_prot[origin_seq][strand] = {}
            strand_cds = origin_cds[strand]
            for i in range(len(strand_cds['order'])):
                cds_id = strand_cds['order'][i]
                seq = strand_cds['details'][cds_id]['seq']
                if str(seq).endswith("TAG"):
                    start = strand_cds['details'][cds_id]['start']
                    end = strand_cds['details'][cds_id]['end']
                    tag_ending_prot[origin_seq][strand][cds_id] = {
                        "start": start,
                        "end": end,
                        "seq": seq,
                        "order_id": i}
                    tag_ending_cds_info[cds_id] = {
                        "start": start,
                        "end": end,
                        "strand": strand,
                        "origin_seq": origin_seq}

    export_csv(
        tag_ending_cds_info,
        tag_ending_cds_info_path,
        ["start", "end", "strand", "origin_seq"])
    return tag_ending_prot, len(tag_ending_cds_info.keys())


def extract_origin_seq(genome_filepath):
    """
    Extract the sequence of the genome
    """
    origin_seqs = {}
    for record in SeqIO.parse(genome_filepath,"fasta"):
        origin_seqs[record.id] = {
            "genome": record,
            "length": len(record.seq),
            "rev_comp_genome": record.reverse_complement()}
    return origin_seqs


def extend_to_next_stop_codon(current_end, genome, next_cds_end):
    """
    Extend a sequence to the next STOP codon on the genome
    """
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
            new_ends.append(new_end)
            if codon != 'TAG':
                to_continue = False
    return new_ends


def extract_possible_seq(start, end, seq, new_ends, strand, origin_seq, 
genome, genome_size):
    """
    Extract the start and ends of possible sequences for a CDS identified as
    potential PYL cds
    """
    pot_pyl_prot_def = {
        "strand": strand,
        "origin_seq": origin_seq,
        "potential_seq": []}

    if strand == "reverse":
        pot_pyl_prot_def["potential_seq"].append({
            "start": genome_size - end + 1,
            "end": genome_size - start + 1,
            "seq": seq})
    else:
        pot_pyl_prot_def["potential_seq"].append({
            "start": start,
            "end": end,
            "seq": seq})
    for new_end in new_ends:
        new_start = start
        new_seq = genome.seq[(start-1):new_end]
        if strand == "reverse":
            new_start = genome_size - new_end + 1
            new_end = genome_size - start + 1
        pot_pyl_prot_def["potential_seq"].append({
            "start": new_start,
            "end": new_end,
            "seq": new_seq})
    return pot_pyl_prot_def


def extract_potential_pyl_proteins(tag_ending_prot, pred_cds, genome_filepath):
    """
    Extract potential PYL proteins
    """
    origin_seqs = extract_origin_seq(genome_filepath)
    pot_pyl_prot_nb = 0
    pot_pyl_prot = {}
    for origin_seq in tag_ending_prot:
        origin_tag_cds = tag_ending_prot[origin_seq]
        origin_pred_cds = pred_cds[origin_seq]
        full_origin_seq = origin_seqs[origin_seq]
        for strand in origin_tag_cds:
            strand_tag_cds = origin_tag_cds[strand]
            strand_pred_cds = origin_pred_cds[strand]
            pred_cds_nb = len(origin_pred_cds[strand]['order'])
            for cds in strand_tag_cds:
                tag_cds = strand_tag_cds[cds]
                start = tag_cds['start']
                end = tag_cds['end']
                seq = tag_cds['seq']
                order_id = tag_cds['order_id']
                genome_size = full_origin_seq["length"]
                if strand == "reverse":
                    previous_start = start
                    start = (genome_size-end)+1
                    end = (genome_size-previous_start+1)
                    genome = full_origin_seq["rev_comp_genome"]
                    next_id = order_id - 1
                    if next_id >= 0:
                        next_cds_id = strand_pred_cds['order'][next_id]
                        next_cds_end = strand_pred_cds['details'][next_cds_id]['start']
                    else:
                        next_cds_end = genome_size
                    next_cds_end = (genome_size-next_cds_end+1)
                else:
                    genome = full_origin_seq["genome"]
                    next_id = order_id + 1
                    if next_id < pred_cds_nb:
                        next_cds_id = strand_pred_cds['order'][next_id]
                        next_cds_end = strand_pred_cds['details'][next_cds_id]['end']
                    else:
                        next_cds_end = genome_size
                new_ends = extend_to_next_stop_codon(end, genome, next_cds_end)
                if len(new_ends) > 0:
                    pot_pyl_prot[cds] = extract_possible_seq(
                        start,
                        end,
                        seq,
                        new_ends,
                        strand,
                        origin_seq,
                        genome,
                        full_origin_seq["length"])
    return pot_pyl_prot


def find_stop_codon_pos_in_seq(string):
    """
    Find position of STOP codon in a sequence
    """
    stop_codon_pos = []
    for i in range(len(string)-1):
        if string.startswith('*', i):
            stop_codon_pos.append(i)
    return stop_codon_pos


def translate(seq):
    """
    Translate a sequence into amino acids
    """
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


def save_potential_pyl_proteins(pot_pyl_prot, pyl_protein_filepath, log_file, 
potential_pyl_seq_info):
    """
    Save potential PYL proteins in a fasta file
    """
    sequences = []
    info = {}
    for prot_id in pot_pyl_prot:
        count = 0
        log_file.write("\t%s\n" % (prot_id))
        info[prot_id] = {
            "start": pot_pyl_prot[prot_id]["potential_seq"][0]["start"],
            "end": pot_pyl_prot[prot_id]["potential_seq"][0]["end"],
            "strand": pot_pyl_prot[prot_id]["strand"],
            "origin_seq": pot_pyl_prot[prot_id]["origin_seq"],
            "alternative_start": "",
            "alternative_end": ""}
        for potential_seq in pot_pyl_prot[prot_id]["potential_seq"]:
            count += 1
            seq_id = "%s_%s" % (prot_id, count)
            desc = "# origin_seq: %s" % (pot_pyl_prot[prot_id]["origin_seq"])
            desc += " # strand: %s" % (pot_pyl_prot[prot_id]["strand"])
            desc += " # start: %s" % (potential_seq["start"])
            desc += " # end: %s" % (potential_seq["end"])
            translated_seq = translate(potential_seq["seq"])
            sequences.append(SeqRecord(
                translated_seq,
                id=seq_id,
                description=desc))
            log_file.write("\t\t%s\t" % (pot_pyl_prot[prot_id]["strand"]))
            log_file.write("%s\t" % (potential_seq["start"]))
            log_file.write("%s\n" % (potential_seq["end"]))
            if count > 1:
                if count > 2:
                    info[prot_id]["alternative_start"] = ",".join([
                        info[prot_id]["alternative_start"],
                        str(potential_seq["start"])])
                    info[prot_id]["alternative_end"] = ",".join([
                        info[prot_id]["alternative_end"],
                        str(potential_seq["end"])])
                else:
                    info[prot_id]["alternative_start"] = str(
                        potential_seq["start"])
                    info[prot_id]["alternative_end"] = str(
                        potential_seq["end"])

    with open(pyl_protein_filepath, "w") as output_file:
        SeqIO.write(sequences, output_file, "fasta")
    export_csv(
        info,
        potential_pyl_seq_info,
        ["start", "end", "strand", "origin_seq", "alternative_start", "alternative_end"])


def predict_pyl_proteins(genome, predicted_cds, pot_pyl_seq, log,
predicted_cds_info, tag_ending_cds_info, potential_pyl_seq_info):
    """
    """
    with open(log, 'w') as log_file:
        pred_cds, pred_cds_nb = extract_predicted_cds(
            predicted_cds,
            predicted_cds_info)
        msg = "Number of predicted CDS: %s\n"  % (pred_cds_nb)
        log_file.write(msg)

        tag_ending_prot, tag_ending_prot_nb = identify_tag_ending_proteins(
            pred_cds,
            tag_ending_cds_info)
        msg = "Number of TAG-ending predicted CDS: %s\n" % (tag_ending_prot_nb)
        log_file.write(msg)

        pot_pyl_prot = extract_potential_pyl_proteins(
            tag_ending_prot,
            pred_cds,
            genome)
        pot_pyl_prot_nb = len(pot_pyl_prot.keys())
        msg = "Number of potential Pyl proteins: %s\n"  % (pot_pyl_prot_nb)
        log_file.write(msg)
        save_potential_pyl_proteins(
            pot_pyl_prot,
            pot_pyl_seq,
            log_file,
            potential_pyl_seq_info)


if __name__ == '__main__':
    predict_pyl_proteins(
        genome=str(snakemake.input.genome),
        predicted_cds=str(snakemake.input.predicted_cds),
        pot_pyl_seq=str(snakemake.output.potential_pyl_sequences),
        log=str(snakemake.output.pyl_protein_prediction_log),
        predicted_cds_info=str(snakemake.output.predicted_cds_info),
        tag_ending_cds_info=str(snakemake.output.tag_ending_cds_info),
        potential_pyl_seq_info=str(snakemake.output.potential_pyl_protein_info))
