import logging

from Bio import SeqIO

from pylprotpredictor import cds
from pylprotpredictor import export


def extract_seqs(seq_filepath):
    """Extract the sequences in a fasta file

    :param seq_filepath: path to a fasta file

    :return: a dictionary with all sequences indexed by their id, their length and their complement sequence
    """
    seqs = {}
    for record in SeqIO.parse(seq_filepath, "fasta"):
        seqs[record.id] = record
    return seqs


def extract_predicted_cds(
        pred_cds_path, pred_cds_info_path, tag_ending_cds_info_path,
        genome_filepath):
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

    for record in SeqIO.parse(pred_cds_path, "fasta"):
        cds_obj = cds.CDS()
        cds_obj.init_from_record(record)
        cds_id = cds_obj.get_id()
        origin_seq_id = cds_obj.get_origin_seq_id()
        strand = cds_obj.get_strand()
        start = cds_obj.get_start()
        end = cds_obj.get_end()

        pred_cds[cds_id] = cds_obj

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

    logging.info("Number of predicted CDS: %s\n" % (len(pred_cds.keys())))
    logging.info("Number of TAG-ending predicted CDS: %s\n" % (len(tag_ending_cds)))

    return pred_cds, ordered_pred_cds, tag_ending_cds


def extract_potential_pyl_cds(
        tag_ending_cds, pred_cds, ordered_pred_cds, pot_pyl_cds_filepath,
        pot_pyl_cds_info_filepath, pot_pyl_cds_obj_filepath):
    """Extract potential PYL CDS from TAG-ending CDS

    :param tag_ending_cds: a list of the ids of the TAG-ending CDS
    :param pred_cds: a dictionary with the predicted CDS represented as CDS objects
    :param ordered_pred_cds: a dictionary with per origin seq, per strand the list of the predicted CDS represented by CDS object
    :param pot_pyl_cds_filepath: path to fasta file in which the protein sequences of the potential PYL CDS are saved
    :param pot_pyl_cds_info_filepath: path to a cvs file to get information about potential PYL CDS
    :param pot_pyl_cds_obj_filepath: path to a JSON file to store the list of potential PYL CDS object
    """
    pot_pyl_cds_nb = 0
    pot_pyl_cds = {}
    info = {}
    sequences = []
    for cds_id in tag_ending_cds:
        cds_obj = pred_cds[cds_id]

        origin_seq_id = cds_obj.get_origin_seq_id()
        strand = cds_obj.get_strand()

        cds_obj.find_next_cds_limit(ordered_pred_cds[origin_seq_id][strand], pred_cds)
        cds_obj.find_alternative_ends()
        cds_obj.extract_possible_alternative_seq()

        if cds_obj.is_potential_pyl_cds():
            info[cds_id] = {
                "start": cds_obj.get_start(),
                "end": cds_obj.get_end(),
                "strand": cds_obj.get_strand(),
                "origin_seq": cds_obj.get_origin_seq_id(),
                "alternative_start": ";".join([str(s) for s in cds_obj.get_alternative_start()]),
                "alternative_end": ";".join([str(s) for s in cds_obj.get_alternative_end()])}
            sequences.append(cds_obj.get_translated_seq())
            sequences += cds_obj.get_translated_alternative_seq()
            pot_pyl_cds.update(cds_obj.export_to_dict())

    logging.info("Number of potential Pyl proteins: %s\n" % (len(pot_pyl_cds.keys())))
    export.export_fasta(sequences, pot_pyl_cds_filepath)
    export.export_csv(
        info,
        pot_pyl_cds_info_filepath,
        ["start", "end", "strand", "origin_seq", "alternative_start", "alternative_end"])
    export.export_json(pot_pyl_cds, pot_pyl_cds_obj_filepath)


def predict_pyl_proteins(
        genome_filepath, pred_cds_filepath, pot_pyl_seq_filepath, log_filepath,
        pred_cds_info_filepath, tag_ending_cds_info_filepath,
        pot_pyl_seq_info_filepath, pot_pyl_cds_obj_filepath):
    """
    Run prediction of potentila PYL CDS:

    - Extraction of predicted CDS into a dictionary
    - Identification of TAG-ending proteins
    - Extraction of potential PYL sequences

    :param genome_filepath: path to file with genome sequence
    :param pred_cds_filepath: path to the output of CDS prediction (Prodigal)
    :param pot_pyl_seq_filepath: path to fasta file with potential PYL CDS sequence
    :param log_filepath: path to log file
    :param pred_cds_info_filepath: path to CSV file with predicted CDS info
    :param tag_ending_cds_info_filepath: path to CSV file with TAG-ending CDS info
    :param pot_pyl_seq_info_filepath: path to CSV file with potential PYL CDS info
    :param pot_pyl_cds_obj_filepath: path to a JSON file to store the list of potential PYL CDS object
    """
    logging.basicConfig(filename=log_filepath, level=logging.INFO, filemode="w")
    pred_cds, ordered_pred_cds, tag_ending_cds = extract_predicted_cds(
        pred_cds_filepath,
        pred_cds_info_filepath,
        tag_ending_cds_info_filepath,
        genome_filepath)
    extract_potential_pyl_cds(
        tag_ending_cds,
        pred_cds,
        ordered_pred_cds,
        pot_pyl_seq_filepath,
        pot_pyl_seq_info_filepath,
        pot_pyl_cds_obj_filepath)
