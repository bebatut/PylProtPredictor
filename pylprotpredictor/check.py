import pandas as pd

from pylprotpredictor import alignment
from pylprotpredictor import cds
from pylprotpredictor import export


def import_cds(pot_pyl_cds_filepath):
    """
    :param pot_pyl_cds_filepath: path to JSON file with collection of potential PYL CDS

    :return: dictionary of the potential PYL CDS
    """
    pot_pyl_cds = {}
    d = export.import_json(pot_pyl_cds_filepath)
    for cds_id in d:
        cds_obj = cds.CDS()
        cds_obj.init_from_dict(d[cds_id])
        pot_pyl_cds[cds_id] = cds_obj
    return pot_pyl_cds


def get_cds_obj(cds_id, pot_pyl_cds):
    """Find the CDS object given an id

    :param cds_id: id of the CDS to find
    :param pot_pyl_cds: dictionary of the potential PYL CDS

    :return: a CDS object
    """
    if cds_id not in pot_pyl_cds:
        raise ValueError("CDS not found for %s" % (cds_id))
    return pot_pyl_cds[cds_id]


def parse_similarity_search_report(pot_pyl_similarity_search, pot_pyl_cds):
    """Parse the similarity search report and add information to the list of
    potential PYL CDS

    :param pot_pyl_similarity_search: path to similarity search report of potential PYL CDS against a reference database
    :param pot_pyl_cds: dictionary of the potential PYL CDS
    """
    similarity_search_report = pd.read_table(
        pot_pyl_similarity_search,
        index_col=None,
        header=None)

    for row in similarity_search_report.itertuples():
        cds_id = row[1].split("-")[0]
        cds_obj = get_cds_obj(cds_id, pot_pyl_cds)

        al = alignment.Alignment()
        al.init_from_search_report_row(row)
        cds_obj.add_id_alignment(row[1], al)

    return pot_pyl_cds


def extract_correct_cds(
        pot_pyl_cds, cons_pot_pyl_seq_filepath, info_filepath):
    """Identify and extract the correct CDS sequence

    :param pot_pyl_cds: dictionary of the potential PYL CDS
    :param cons_pot_pyl_seq_filepath: path to a FASTA file for the conserved CDS
    :param info_filepath: path to a CSV file with final information about the CDS
    """
    cons_sequences = []
    info = {}
    keys = list(pot_pyl_cds.keys())
    keys.sort()
    for cds_id in keys:
        cds_obj = pot_pyl_cds[cds_id]
        cds_obj.identify_cons_rej_cds()

        cons_seq = cds_obj.get_conserved_cds()
        cons_sequences.append(cons_seq.get_seqrecord())

        rej_seq = cds_obj.get_rejected_cds()
        info[cds_id] = {
            "conserved_start": cons_seq.get_start(),
            "conserved_end": cons_seq.get_end(),
            "strand": cds_obj.get_strand(),
            "origin_seq": cds_obj.get_origin_seq_id(),
            "original_start": cds_obj.get_start(),
            "original_end": cds_obj.get_end(),
            "rejected_start": ";".join([str(s.get_start()) for s in rej_seq]),
            "rejected_end": ";".join([str(s.get_end()) for s in rej_seq])}

    export.export_fasta(cons_sequences, cons_pot_pyl_seq_filepath)
    export.export_csv(
        info,
        info_filepath,
        ["conserved_start", "conserved_end", "strand", "origin_seq", "original_start", "original_end", "rejected_start", "rejected_end"])


def check_pyl_proteins(
        pot_pyl_similarity_search, pot_pyl_cds_filepath,
        cons_pot_pyl_seq, info_filepath):
    """
    Check predicted PYL CDS:

    - Get the potential PYL CDS
    - Parse the similarity search report
    - Identify and extract the correct CDS sequence (the one with the lowest evalue)

    :param pot_pyl_similarity_search: path to similarity search report of potential PYL CDS against a reference database
    :param pot_pyl_cds_filepath: path to JSON file with collection of potential PYL CDS
    :param cons_pot_pyl_seq_filepath: path to a FASTA file for the conserved CDS
    :param info_filepath: path to a CSV file with final information about the CDS
    """
    pot_pyl_cds = import_cds(pot_pyl_cds_filepath)
    parse_similarity_search_report(pot_pyl_similarity_search, pot_pyl_cds)
    extract_correct_cds(pot_pyl_cds, cons_pot_pyl_seq, info_filepath)
