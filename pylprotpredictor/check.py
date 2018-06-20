import argparse
import pandas as pd

from pathlib import Path

try:
    from pylprotpredictor import alignment
    from pylprotpredictor import cds
    from pylprotpredictor import export
except ImportError:
    from . import alignment
    from . import cds
    from . import export


def import_cds(cds_obj_filepath):
    """
    :param cds_obj_filepath: path to JSON file with collection of CDS objects

    :return: dictionary of the CDS objects
    """
    cds_objects = {}
    d = export.import_json(cds_obj_filepath)
    for cds_id in d:
        cds_obj = cds.CDS()
        cds_obj.init_from_dict(d[cds_id])
        cds_objects[cds_id] = cds_obj
    return cds_objects


def get_cds_obj(cds_id, pred_cds):
    """Find the CDS object given an id

    :param cds_id: id of the CDS to find
    :param pred_cds: dictionary of the predicted CDS

    :return: a CDS object
    """
    if cds_id not in pred_cds:
        raise ValueError("CDS not found for %s" % (cds_id))
    return pred_cds[cds_id]


def parse_similarity_search_report(pot_pyl_similarity_search, pred_cds):
    """Parse the similarity search report and add information to the list of
    potential PYL CDS

    :param pot_pyl_similarity_search: path to similarity search report of potential PYL CDS against a reference database
    :param pred_cds: dictionary of the predicted CDS
    """
    similarity_search_report = pd.read_table(
        pot_pyl_similarity_search,
        index_col=None,
        header=None)

    for row in similarity_search_report.itertuples():
        cds_id = row[1].split("-")[0]
        cds_obj = get_cds_obj(cds_id, pred_cds)

        al = alignment.Alignment()
        al.init_from_search_report_row(row)
        cds_obj.add_id_alignment(row[1], al)

    return pred_cds


def extract_correct_cds(pred_cds, cons_pred_cds_seq, info_filepath):
    """Identify and extract the correct CDS sequence

    :param pred_cds: dictionary of the predicted CDS
    :param cons_pred_cds_seq: path to a FASTA file for the conserved CDS sequences
    :param info_filepath: path to a CSV file with final information about the CDS
    """
    cons_cds_sequences = []
    info = {}
    keys = list(pred_cds.keys())
    keys.sort()

    for cds_id in keys:
        cds_obj = pred_cds[cds_id]
        cons_al = alignment.Alignment()

        if cds_obj.is_potential_pyl():
            cons_al = cds_obj.identify_cons_rej_cds()

        cons_seq = cds_obj.get_conserved_cds()
        rej_seq = cds_obj.get_rejected_cds()

        if cons_seq is None:
            cons_seq = cds_obj

        cons_cds_sequences.append(cons_seq.get_seqrecord())

        info[cds_id] = {
            "status": cds_obj.get_status(),
            "conserved_start": cons_seq.get_start(),
            "conserved_end": cons_seq.get_end(),
            "strand": cds_obj.get_strand(),
            "conserved_stop_codon": cons_seq.get_stop_codon(),
            "origin_seq": cds_obj.get_origin_seq_id(),
            "original_start": cds_obj.get_start(),
            "original_end": cds_obj.get_end(),
            "original_stop_codon": cds_obj.get_stop_codon(),
            "rejected_start": ";".join([str(s.get_start()) for s in rej_seq]),
            "rejected_end": ";".join([str(s.get_end()) for s in rej_seq]),
            "matched_RefSeq90": cons_al.get_sseqid()}

    export.export_fasta(cons_cds_sequences, cons_pred_cds_seq)
    export.export_csv(
        info,
        info_filepath,
        ["status", "conserved_start", "conserved_end", "strand", "conserved_stop_codon", "origin_seq", "original_start", "original_end", "original_stop_codon", "rejected_start", "rejected_end", "matched_RefSeq90"])


def check_pyl_proteins(
        pot_pyl_similarity_search, pred_cds_obj_filepath,
        cons_pred_cds_seq, info_filepath):
    """
    Check predicted PYL CDS:

    - Get the potential PYL CDS
    - Parse the similarity search report
    - Identify and extract the correct CDS sequence (the one with the lowest evalue and longest alignment for potential PYL)

    :param pot_pyl_similarity_search: path to similarity search report of potential PYL CDS against a reference database
    :param pred_cds_obj_filepath: path to generated JSON file to store the list of predicted CDS objects
    :param cons_pred_cds_seq: path to a FASTA file for the conserved CDS sequences
    :param info_filepath: path to a CSV file with final information about the CDS
    """
    pred_cds = import_cds(pred_cds_obj_filepath)
    parse_similarity_search_report(pot_pyl_similarity_search, pred_cds)
    extract_correct_cds(pred_cds, cons_pred_cds_seq, info_filepath)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Check predicted PYL CDS')
    parser.add_argument('--pot_pyl_similarity_search', help='path to similarity search report of potential PYL CDS against a reference database')
    parser.add_argument('--pred_cds_obj_filepath', help='path to generated JSON file to store the list of predicted CDS objects')
    parser.add_argument('--cons_pred_cds_seq_filepath', help='path to a FASTA file for the conserved CDS sequences')
    parser.add_argument('--info_filepath', help='path to a CSV file with final information about the CDS')

    args = parser.parse_args()

    check_pyl_proteins(
        pot_pyl_similarity_search=Path(args.pot_pyl_similarity_search),
        pred_cds_obj_filepath=Path(args.pred_cds_obj_filepath),
        cons_pred_cds_seq=Path(args.cons_pred_cds_seq_filepath),
        info_filepath=Path(args.info_filepath))
