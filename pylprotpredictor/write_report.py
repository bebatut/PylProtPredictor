import pandas as pd
from snakemake.utils import report


def extract_row_number(csv_filepath):
    """
    Extract row number of a CSV file
    """
    df = pd.read_csv(csv_filepath, index_col=0, header=0)
    return df.shape[0]


if __name__ == '__main__':
    predicted_cds_nb = extract_row_number(
        str(snakemake.input.predicted_cds_info))
    tag_ending_cds_nb = extract_row_number(
        str(snakemake.input.tag_ending_cds))
    potential_pyl_protein_nb = extract_row_number(
        str(snakemake.input.potential_pyl_proteins))
    cons_pot_pyl_seq_nb = extract_row_number(
        str(snakemake.input.conserved_potential_pyl_sequences))

    report("""
    Prediction of PYL proteins
    ==========================

    {predicted_cds_nb} sequences are predicted as CDS (see Table
    predicted_cds_info_).

    From these {predicted_cds_nb} prediced CDS, {tag_ending_cds_nb} are ending
    with a TAG STOP codon (see Table tag_ending_cds_).
    {potential_pyl_protein_nb} are predicted as potential PYL sequences (see
    Table potential_pyl_proteins_).

    After similarity search, {cons_pot_pyl_seq_nb} sequences (of
    {potential_pyl_protein_nb} predicted PYL sequences) are conserved as
    potential PYL proteins (see Table conserved_potential_pyl_sequences_).
    """, str(snakemake.output.report), **snakemake.input)
