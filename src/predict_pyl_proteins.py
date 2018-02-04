#!/usr/bin/env python

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd
from src import export
from src import cds
import logging


def get_string_seq(seq):
    """Return the string of the sequence

    :param seq: SeqRecord object

    :return: string corresponding to the sequence
    """
    return str(seq.seq)


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


def extract_seqs(seq_filepath):
    """Extract the sequences in a fasta file

    :param seq_filepath: path to a fasta file

    :return: a dictionary with all sequences indexed by their id, their length and their complement sequence
    """
    seqs = {}
    for record in SeqIO.parse(seq_filepath,"fasta"):
        seqs[record.id] = record
    return seqs


def extract_predicted_cds(pred_cds_path, pred_cds_info_path, tag_ending_cds_info_path, genome_filepath):
    """Extract the list of predicted CDS and identify the CDS ending with TAG STOP codon

    :param pred_cds_path: path to the output of CDS prediction (Prodigal)
    :param pred_cds_info_path: path to a CSV file in which the information (start, end, strand, origin) are collected for each predicted CDS
    :param tag_ending_cds_info_path: path to CSV file to export the information about the TAG ending CDS
    :param genome_filepath: path to reference genome

    :return: a dictionary with the predicted CDS represented by CDS object
    :return: a dictionary with per origin seq, per strand the list of the predicted CDS represented by CDS object
    :return: a list of the id of TAG-ending CDS
    """
    pred_cds = {}
    ordered_pred_cds = {}
    pred_cds_info = {}
    tag_ending_cds = []
    tag_ending_cds_info = {}

    origin_seqs = extract_seqs(genome_filepath)

    for record in SeqIO.parse(pred_cds_path,"fasta"):
        cds_obj = cds.CDS()
        cds_obj.init_from_record(record)
        cds_id = cds_obj.get_id()
        origin_seq_id = cds_obj.get_origin_seq_id()
        strand = cds_obj.get_strand()
        start = cds_obj.get_start()
        end = cds_obj.get_end()
        # 
        pred_cds[cds_id] = cds_obj
        # 
        ordered_pred_cds.setdefault(origin_seq_id, {})
        ordered_pred_cds[origin_seq_id].setdefault(strand, [])

        cds_obj.set_order(len(ordered_pred_cds[origin_seq_id][strand]))
        ordered_pred_cds[origin_seq_id][strand].append(cds_id)
        
        pred_cds_info[cds_id] = {"start": start, "end": end, "strand": strand, "origin_seq": origin_seq_id}

        if cds_obj.is_tag_ending_seq():
            tag_ending_cds.append(cds_id)
            tag_ending_cds_info[cds_id] = {"start": start, "end": end, "strand": strand, "origin_seq": origin_seq_id}
            cds_obj.set_origin_seq(origin_seqs[origin_seq_id])
    
    export.export_csv(pred_cds_info, pred_cds_info_path, ["start", "end", "strand", "origin_seq"])
    export.export_csv(tag_ending_cds_info, tag_ending_cds_info_path, ["start", "end", "strand", "origin_seq"])

    logging.info("Number of predicted CDS: %s\n"  % (len(pred_cds.keys())))
    logging.info("Number of TAG-ending predicted CDS: %s\n" % (len(tag_ending_cds)))
    
    return pred_cds, ordered_pred_cds, tag_ending_cds


def extract_potential_pyl_cds(tag_ending_cds, pred_cds, ordered_pred_cds):
    """Extract potential PYL CDS from TAG-ending CDS

    :param tag_ending_cds: a list of the ids of the TAG-ending CDS
    :param pred_cds: a dictionary with the predicted CDS represented as CDS objects

    :return: dictionary with for each potential PYL CDS (keys) a dictionary with several information and the description of the possible sequences for the CDS
    """
    pot_pyl_cds = {}
    for cds_id in tag_ending_cds:
        cds_obj = pred_cds[cds_id]

        origin_seq_id = cds_obj.get_origin()
        strand = cds_obj.get_strand()

        cds_obj.find_next_cds_limit(ordered_pred_cds[origin_seq_id][strand], pred_cds)
        cds_obj.find_alternative_ends()
        cds_obj.extract_possible_seq()

    logging.info("Number of potential Pyl proteins: %s\n"  % (len(pot_pyl_prot.keys())))
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
            seq_id = "%s_%s" % (cds_id, count)
            desc = "# origin_seq: %s # strand: %s # start: %s # end: %s" % (
                pot_pyl_cds[cds_id]["origin_seq"],
                pot_pyl_cds[cds_id]["strand"],
                potential_seq["start"],
                potential_seq["end"])
            translated_seq = translate(Seq(potential_seq["seq"]))
            sequences.append(SeqRecord(translated_seq, id=cds_id, description=desc))
            logging.info("\t\t%s\t" % (pot_pyl_cds[cds_id]["strand"]))
            logging.info("%s\t" % (potential_seq["start"]))
            logging.info("%s\n" % (potential_seq["end"]))
            
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

    export.export_fasta(sequences, pot_pyl_cds_filepath)
    export.export_csv(
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
    logging.basicConfig(filename=log, level=logging.INFO, filemode="w")
    pred_cds, pred_cds_nb, tag_ending_cds, tag_ending_cds_nb = extract_predicted_cds(predicted_cds, predicted_cds_info, log_file)
    pot_pyl_prot = extract_potential_pyl_cds(tag_ending_cds, pred_cds, genome)
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
