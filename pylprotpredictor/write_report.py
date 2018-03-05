import argparse
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

    :param pred_cds: path to CSV file with predicted CDS info
    :param tag_ending_cds: path to CSV file with TAG-ending CDS info
    :param pot_pyl_cds: path to CSV file with potential PYL CDS info
    :param final_cds_info: path to a CSV file with final information about the CDS
    :param report_filepath: path to HTML file in which writing the report
    """
    pred_cds_nb = extract_row_number(pred_cds)
    tag_ending_cds_nb = extract_row_number(tag_ending_cds)
    pot_pyl_cds_nb = extract_row_number(pot_pyl_cds)

    df = pd.read_csv(final_cds, index_col=0, header=0)
    rej_nb = len(df.query("conserved_start == original_start and conserved_end == original_end"))
    final_cds_nb = pot_pyl_cds_nb - rej_nb

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
        {pot_pyl_cds_nb} predicted PYL sequences) are conserved ({rej_nb}
        rejected) as potential PYL proteins (see Table final_cds_).
        """,
        report_filepath,
        predicted_cds=pred_cds,
        tag_ending_cds=tag_ending_cds,
        potential_pyl_cds=pot_pyl_cds,
        final_cds=final_cds)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Write HTML report to summarize the full analysis')
    parser.add_argument('--pred_cds', help='path to CSV file with predicted CDS info')
    parser.add_argument('--tag_ending_cds', help='path to CSV file with TAG-ending CDS info')
    parser.add_argument('--pot_pyl_cds', help='path to CSV file with potential PYL CDS info')
    parser.add_argument('--final_cds_info', help='path to a CSV file with final information about the CDS')
    parser.add_argument('--report_filepath', help='path to HTML file in which writing the report')

    args = parser.parse_args()

    write_report(
        pred_cds=args.pred_cds,
        tag_ending_cds=args.tag_ending_cds,
        pot_pyl_cds=args.pot_pyl_cds,
        final_cds=args.final_cds_info,
        report_filepath=args.report_filepath)
