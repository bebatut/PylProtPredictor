import pandas as pd

from snakemake.utils import report


def extract_row_number(csv_filepath):
    """Extract row number of a CSV file

    :param csv_filepath: path to a CSV file

    :return: an integer corresponding to the number of lines in the CSV file
    """
    df = pd.read_csv(csv_filepath, index_col=0, header=0)
    return df.shape[0]


def write_report(pred_cds, tag_ending_cds, pot_pyl_cds, final_cds, report_filepath):
    """Write HTML report to summarize the full analysis

    :param pred_cds:
    :param tag_ending_cds:
    :param pot_pyl_cds:
    :param final_cds_info:
    :param report_filepath:
    """
    pred_cds_nb = extract_row_number(pred_cds)
    tag_ending_cds_nb = extract_row_number(tag_ending_cds)
    pot_pyl_cds_nb = extract_row_number(pot_pyl_cds)
    final_cds_nb = extract_row_number(final_cds)

    print("%s %s %s %s" % (pred_cds_nb, tag_ending_cds_nb, pot_pyl_cds_nb, final_cds_nb))

    report(
        """
        Prediction of PYL proteins
        ==========================

        {pred_cds_nb} sequences are predicted as CDS (see Table
        predicted_cds_).

        From these {pred_cds_nb} prediced CDS, {tag_ending_cds_nb} are ending
        with a TAG STOP codon (see Table tag_ending_cds_).
        {pot_pyl_cds_nb} are predicted as potential PYL sequences (see
        Table potential_pyl_cds_).

        After similarity search, {final_cds_nb} sequences (of
        {pot_pyl_cds_nb} predicted PYL sequences) are conserved as
        potential PYL proteins (see Table final_cds_).
        """,
        report_filepath,
        predicted_cds=pred_cds,
        tag_ending_cds=tag_ending_cds,
        potential_pyl_cds=pot_pyl_cds,
        final_cds=final_cds)
