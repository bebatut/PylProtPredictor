#!/usr/bin/env python

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Data import CodonTable
import pandas as pd


def transform_strand(strand_id):
    """Transform strand from numerical value to string value

    :param strand_id: numerical value to represent a strand (1 or -1)

    :return: string value (forward or reverse) for the strand
    """
    if strand_id == "-1":
        return 'reverse'
    elif strand_id == "1":
        return 'forward'
    else:
        raise ValueError("Wrong strand_id: {}".format(strand_id))


def get_string_seq(seq):
    """Return the string of the sequence

    :param seq: SeqRecord object

    :return: string corresponding to the sequence
    """
    return str(seq.seq)


def extract_seq_desc(desc):
    """Extract from description the seq id, the origin sequence, start, end and
    strand from a predicted CDS

    :param desc: description of a prediced CDS with Prodigal

    :return: id of predicted CDS
    :return: id of the origin sequence
    :return: start position of the predicted CDS
    :return: end position of the predicted CDS
    :return: strand of the predicted CDS
    """
    split_description = desc.split("#")
    seq_id = split_description[0][:-1]
    origin_seq = "_".join(seq_id.split("_")[:-1])
    if seq_id.find("|") != -1:
        seq_id = "cds_%s" % (seq_id.split("_")[-1])
    start = int(split_description[1].replace(" ",""))
    end = int(split_description[2].replace(" ",""))
    strand = transform_strand(split_description[3].replace(" ",""))
    return seq_id, origin_seq, start, end, strand


def translate(seq):
    """Translate a sequence into amino acids while replacing any possible STOP 
    codon encoded by TAG by a Pyl amino acid

    :param seq: a SeqRecord object

    :return: string with the corresponding amino acid sequence with the TAG encoded STOP are replaced by Pyl amino acid
    """
    translated_seq = seq.translate()
    str_seq = str(seq)
    n = 3
    codons = [str_seq[i:i+n] for i in range(0, len(str_seq), n)]
    for i in find_stop_codon_pos_in_seq(str(translated_seq)):
        if codons[i] != "TAG":
            raise ValueError("Stop codon (different from TAG) found inside the sequence")
        mutable_seq = translated_seq.tomutable()
        mutable_seq[i] = "O"
        translated_seq = mutable_seq.toseq()
    return translated_seq


def find_stop_codon_pos_in_seq(seq):
    """Find position of STOP codon in a sequence

    :param seq: string sequence of amino acids

    :return: list of position for possible STOP codons in a sequence
    """
    stop_codon_pos = []
    for i in range(len(seq)-1):
        if seq.startswith('*', i):
            stop_codon_pos.append(i)
    return stop_codon_pos


def export_csv(dictionary, csv_filepath, col):
    """Export a dictionary into a csv file with specific columns

    :param dictionary: dictionary with sequence ids as keys and start, end, strand, and origin sequences as values
    :param csv_filepath: path to CSV file to export
    :param col: list of columns (keys of dictionary's subdict) to conserve
    """
    dict_df = pd.DataFrame(data=dictionary).transpose()
    dict_df = dict_df[col]
    dict_df.to_csv(csv_filepath)


def extract_predicted_cds(pred_cds_path, pred_cds_info_path):
    """Extract the list of predicted CDS

    :param pred_cds_path: path to the output of CDS prediction (Prodigal)
    :param pred_cds_info_path: path to a CSV file in which the information (start, end, strand, origin) are collected for each predicted CDS

    :return: a dictionary with per origin seq, per strand are detailed the predicted CDS (start, end, sequence) and their order
    :return: number of predicted CDS
    """
    pred_cds = {}
    pred_cds_info = {}
    for record in SeqIO.parse(pred_cds_path,"fasta"):
        seq_id, origin_seq, start, end, strand = extract_seq_desc(record.description)
        pred_cds.setdefault(origin_seq, {})
        pred_cds[origin_seq].setdefault(strand, {'details': {}, 'order': []})
        pred_cds[origin_seq][strand]['details'][seq_id] = {'start': start, 'end': end, 'seq': record.seq}
        pred_cds[origin_seq][strand]['order'].append(seq_id)
        pred_cds_info[seq_id] = {"start": start, "end": end, "strand": strand, "origin_seq": origin_seq}
    export_csv(pred_cds_info, pred_cds_info_path, ["start", "end", "strand", "origin_seq"])
    return pred_cds, len(pred_cds_info.keys())


def extract_seqs(seq_filepath):
    """Extract the sequences in a fasta file

    :param seq_filepath: path to a fasta file

    :return: a dictionary with all sequences indexed by their id, their length and their complement sequence
    """
    seqs = {}
    for record in SeqIO.parse(seq_filepath,"fasta"):
        seqs[record.id] = {
            "genome": record,
            "rev_comp_genome": record.reverse_complement()}
    return seqs


def identify_tag_ending_cds(pred_cds, tag_ending_cds_info_path):
    """Identify CDS ending with TAG STOP codon

    :param pred_cds: a dictionary with per origin seq, per strand are detailed the predicted CDS (start, end, sequence) and their order
    :param tag_ending_cds_info_path: path to CSV file to export the information about the TAG ending CDS

    :return: dictionary with with per origin seq, per strand are detailed the TAG-ending CDS (start, end, sequence) and their order
    :return: number of TAG-ending CDS
    """
    tag_ending_cds = {}
    tag_ending_cds_info = {}
    for origin_seq in pred_cds:
        origin_cds = pred_cds[origin_seq]
        for strand in origin_cds:
            strand_cds = origin_cds[strand]
            for i in range(len(strand_cds['order'])):
                cds_id = strand_cds['order'][i]
                seq = strand_cds['details'][cds_id]['seq']
                if str(seq).endswith("TAG"):
                    tag_ending_cds.setdefault(origin_seq, {})
                    tag_ending_cds[origin_seq].setdefault(strand, {})
                    start = strand_cds['details'][cds_id]['start']
                    end = strand_cds['details'][cds_id]['end']
                    tag_ending_cds[origin_seq][strand][cds_id] = {
                        "start": start,
                        "end": end,
                        "seq": seq,
                        "order_id": i}
                    tag_ending_cds_info[cds_id] = {
                        "start": start,
                        "end": end,
                        "strand": strand,
                        "origin_seq": origin_seq}
    export_csv(tag_ending_cds_info, tag_ending_cds_info_path, ["start", "end", "strand", "origin_seq"])
    return tag_ending_cds, len(tag_ending_cds_info.keys())


def find_alternative_ends(current_end, genome, next_cds_end):
    """
    Find alternative ends (on the same ORF) for a CDS until the next found STOP 
    codon on the genome (or its complement if the CDS is on the reverse strand)

    :param current_end: position of the current end of the CDS on the genome
    :param genome: sequence string of the reference genome
    :param next_cds_end: position of the end of the next CDS on the genome

    :return: list of possible alternative ends (STOP codons after on the genome but before the end of the next CDS) for the CDS
    """
    genome_size = len(genome)
    stop_codons = CodonTable.unambiguous_dna_by_id[1].stop_codons
    new_end = current_end
    to_continue = (new_end+3) < genome_size and (new_end+3) < next_cds_end 
    new_ends = []
    while to_continue:
        codon = str(genome[new_end:(new_end + 3)])
        if codon not in stop_codons:
            new_end += 3
            to_continue = (new_end+3) < genome_size and (new_end+3) < next_cds_end
        else:
            new_end += 3
            new_ends.append(new_end)
            if codon != 'TAG':
                to_continue = False
            else:
                to_continue = (new_end+3) < genome_size and (new_end+3) < next_cds_end
    return new_ends


def extract_possible_seq(start, end, new_ends, strand, origin_seq_id, origin_seq):
    """
    Extract the start, end and sequence of different possible sequences for a CDS identified as
    potential PYL CDS

    :param start: current start of the CDS
    :param end: current end of the CDS
    :param new_ends: list of possible new ends (next TAG or next codon STOP) for the extended sequences
    :param strand: strand of the CDS
    :param origin_seq_id: id of the origin sequence
    :param origin_seq: sequence string of the origin sequence

    :return: a dictionary with several information and the description of the possible sequences for the CDS
    """
    pot_pyl_cds = {"strand": strand, "origin_seq": origin_seq_id, "potential_seq": []}
    origin_seq_size = len(origin_seq)
    seq = origin_seq[(start-1):end]
    if strand == "reverse":
        pot_pyl_cds["potential_seq"].append({
            "start": origin_seq_size - end + 1,
            "end": origin_seq_size - start + 1,
            "seq": seq})
    else:
        pot_pyl_cds["potential_seq"].append({"start": start, "end": end, "seq": seq})
    for new_end in new_ends:
        new_start = start
        new_seq = origin_seq[(start-1):new_end]
        if strand == "reverse":
            new_start = origin_seq_size - new_end + 1
            new_end = origin_seq_size - start + 1
        pot_pyl_cds["potential_seq"].append({"start": new_start, "end": new_end, "seq": new_seq})
    return pot_pyl_cds


def extract_potential_pyl_cds(tag_ending_cds, pred_cds, genome_filepath):
    """Extract potential PYL CDS from TAG-ending CDS

    :param tag_ending_cds: dictionary with with per origin seq, per strand are detailed the TAG-ending CDS (start, end, sequence) and their order
    :param pred_cds: a dictionary with per origin seq, per strand are detailed the predicted CDS (start, end, sequence) and their order
    :param genome_filepath: path to reference genome

    :return: dictionary with for each potential PYL CDS (keys) a dictionary with several information and the description of the possible sequences for the CDS
    """
    origin_seqs = extract_seqs(genome_filepath)
    pot_pyl_cds = {}
    for origin_seq_id in tag_ending_cds:
        origin_tag_cds = tag_ending_cds[origin_seq_id]
        origin_pred_cds = pred_cds[origin_seq_id]
        full_origin_seq = origin_seqs[origin_seq_id]
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
                origin_seq = get_string_seq(full_origin_seq["genome"])
                origin_seq_size = len(origin_seq)
                if strand == "reverse":
                    previous_start = start
                    start = (origin_seq_size-end)+1
                    end = (origin_seq_size-previous_start+1)
                    origin_seq = get_string_seq(full_origin_seq["rev_comp_genome"])
                    next_id = order_id - 1
                    if next_id >= 0:
                        next_cds_id = strand_pred_cds['order'][next_id]
                        next_cds_end = strand_pred_cds['details'][next_cds_id]['start']
                    else:
                        next_cds_end = 0
                    next_cds_end = (origin_seq_size-next_cds_end+1)
                else:
                    next_id = order_id + 1
                    if next_id < pred_cds_nb:
                        next_cds_id = strand_pred_cds['order'][next_id]
                        next_cds_end = strand_pred_cds['details'][next_cds_id]['end']
                    else:
                        next_cds_end = origin_seq_size
                new_ends = find_alternative_ends(end, origin_seq, next_cds_end)
                if len(new_ends) > 0:
                    pot_pyl_cds[cds] = extract_possible_seq(
                        start,
                        end,
                        new_ends,
                        strand,
                        origin_seq_id,
                        origin_seq)
    return pot_pyl_cds


def save_potential_pyl_cds(pot_pyl_cds, pot_pyl_cds_filepath, log_file, pot_pyl_cds_info_filepath):
    """Save protein sequence of the potential PYL CDS in a fasta file

    :param pot_pyl_cds: dictionary with for each potential PYL CDS (keys) a dictionary with several information and the description of the possible sequences for the CDS
    :param pot_pyl_cds_filepath: path to fasta file in which the protein sequences of the potential PYL CDS are saved
    :param log_file: open stream on a log file
    :param pot_pyl_cds_info_filepath: path to a cvs file to get information about potential PYL CDS
    """
    sequences = []
    info = {}
    
    for cds_id in pot_pyl_cds:
        count = 0
        log_file.write("\t%s\n" % (cds_id))
        info[cds_id] = {
            "start": pot_pyl_cds[cds_id]["potential_seq"][0]["start"],
            "end": pot_pyl_cds[cds_id]["potential_seq"][0]["end"],
            "strand": pot_pyl_cds[cds_id]["strand"],
            "origin_seq": pot_pyl_cds[cds_id]["origin_seq"],
            "alternative_start": "",
            "alternative_end": ""}
        
        for potential_seq in pot_pyl_cds[cds_id]["potential_seq"]:
            count += 1
            desc = "%s_%s # origin_seq: %s # strand: %s # start: %s # end: %s" % (
                cds_id, 
                count,
                pot_pyl_cds[cds_id]["origin_seq"],
                pot_pyl_cds[cds_id]["strand"],
                potential_seq["start"],
                potential_seq["end"])
            translated_seq = translate(Seq(potential_seq["seq"]))
            sequences.append(SeqRecord(translated_seq, id=cds_id, description=desc))
            log_file.write("\t\t%s\t" % (pot_pyl_cds[cds_id]["strand"]))
            log_file.write("%s\t" % (potential_seq["start"]))
            log_file.write("%s\n" % (potential_seq["end"]))
            
            if count == 1:
                info[cds_id]["alternative_start"] = str(potential_seq["start"])
                info[cds_id]["alternative_end"] = str(potential_seq["end"])
            else:
                info[cds_id]["alternative_start"] = ",".join([
                    info[cds_id]["alternative_start"],
                    str(potential_seq["start"])])
                info[cds_id]["alternative_end"] = ",".join([
                    info[cds_id]["alternative_end"],
                    str(potential_seq["end"])])

    with open(pot_pyl_cds_filepath, "w") as output_file:
        SeqIO.write(sequences, output_file, "fasta")
    export_csv(
        info,
        pot_pyl_cds_info_filepath,
        ["start", "end", "strand", "origin_seq", "alternative_start", "alternative_end"])


def predict_pyl_proteins(genome, predicted_cds, pot_pyl_seq, log,
predicted_cds_info, tag_ending_cds_info, potential_pyl_seq_info):
    """
    Run prediction of PYL proteins:

    - Extraction of predicted CDS into a dictionary
    - Identification of TAG-ending proteins
    - Extraction of potential PYL sequences
    
    :param genome:
    :param predicted_cds:
    :param pot_pyl_seq:
    :param predicted_cds_info:
    :param tag_ending_cds_info:
    :param potential_pyl_seq_info:
    """
    with open(log, 'w') as log_file:
        pred_cds, pred_cds_nb = extract_predicted_cds(predicted_cds, predicted_cds_info)
        msg = "Number of predicted CDS: %s\n"  % (pred_cds_nb)
        log_file.write(msg)
        tag_ending_cds, tag_ending_cds_nb = identify_tag_ending_cds(pred_cds, tag_ending_cds_info)
        msg = "Number of TAG-ending predicted CDS: %s\n" % (tag_ending_cds_nb)
        log_file.write(msg)
        pot_pyl_prot = extract_potential_pyl_cds(tag_ending_cds, pred_cds, genome)
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
