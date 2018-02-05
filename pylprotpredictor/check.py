import pandas as pd

from Bio import SeqIO


def parse_similarity_search_report(potential_pyl_similarity_search):
    """
    Parse the similarity search report
    """
    similarity_search_report = pd.read_table(
        potential_pyl_similarity_search,
        index_col=None,
        header=None)
    cds_report = {}
    for row in similarity_search_report.itertuples():
        seq_id = row[1]
        cds_id = "_".join(seq_id.split("_")[:-1])
        evalue = float(row[11])
        cds_report.setdefault(cds_id, {
            "conserved_seq": "",
            "rejected_seq": [],
            "evalue": 10})
        if evalue < cds_report[cds_id]["evalue"]:
            if cds_report[cds_id]["conserved_seq"] != "":
                cds_report[cds_id]["rejected_seq"].append(
                    cds_report[cds_id]["conserved_seq"])
            cds_report[cds_id]["conserved_seq"] = seq_id
            cds_report[cds_id]["evalue"] = evalue
        else:
            cds_report[cds_id]["rejected_seq"].append(seq_id)
    return cds_report


def extract_seq_info(description):
    """
    """
    split_desc = description.split(" # ")
    start = split_desc[3].split(": ")[-1]
    end = split_desc[4].split(": ")[-1]
    strand = split_desc[2].split(": ")[-1]
    return start, end, strand


def extract_conserved_rejected_sequences(cds_report, potential_pyl_seq):
    """
    """
    conserved_seq = []
    cons_seq_info = {}
    rej_seq_info = {}
    for record in SeqIO.parse(potential_pyl_seq, "fasta"):
        seq_id = record.id
        cds_id = "_".join(seq_id.split("_")[:-1])
        seq_nb = seq_id.split("_")[-1]
        start, end, strand = extract_seq_info(record.description)
        if cds_id in cds_report:
            if seq_nb == "1":
                if seq_id in cds_report[cds_id]["conserved_seq"]:
                    rej_seq_info.setdefault(cds_id, {
                        "start": "",
                        "end": "",
                        "strand": "",
                        "rejected_alternative_start": "",
                        "rejected_alternative_end": "",
                        "comment": ""})
                    rej_seq_info[cds_id]["start"] = start
                    rej_seq_info[cds_id]["end"] = end
                    rej_seq_info[cds_id]["strand"] = strand
                else:
                    cons_seq_info.setdefault(cds_id, {
                        "start": "",
                        "end": "",
                        "strand": "",
                        "original_start": "",
                        "original_end": ""})
                    cons_seq_info[cds_id]["original_start"] = start
                    cons_seq_info[cds_id]["original_end"] = end
            else:
                if seq_id in cds_report[cds_id]["conserved_seq"]:
                    cons_seq_info.setdefault(cds_id, {
                        "start": "",
                        "end": "",
                        "strand": "",
                        "original_start": "",
                        "original_end": ""})
                    conserved_seq.append(record)
                    cons_seq_info[cds_id]["start"] = start
                    cons_seq_info[cds_id]["end"] = end
                    cons_seq_info[cds_id]["strand"] = strand
                else:
                    rej_seq_info.setdefault(cds_id, {
                        "start": "",
                        "end": "",
                        "strand": "",
                        "rejected_alternative_start": "",
                        "rejected_alternative_end": "",
                        "comment": ""})
                    if rej_seq_info[cds_id]["rejected_alternative_start"] != "":
                        rej_seq_info[cds_id]["rejected_alternative_start"] = ",".join([
                            rej_seq_info[cds_id]["rejected_alternative_start"],
                            start])
                        rej_seq_info[cds_id]["rejected_alternative_end"] = ",".join([
                            rej_seq_info[cds_id]["rejected_alternative_end"],
                            end])
                    else:
                        rej_seq_info[cds_id]["rejected_alternative_start"] = start
                        rej_seq_info[cds_id]["rejected_alternative_end"] = end
        else:
            rej_seq_info.setdefault(cds_id, {
                "start": "",
                "end": "",
                "strand": "",
                "rejected_alternative_start": "",
                "rejected_alternative_end": "",
                "comment": ""})
            rej_seq_info[cds_id]["comment"] = "None of the possible sequences \
            match the reference database"
            if seq_nb != "1":
                if rej_seq_info[cds_id]["rejected_alternative_start"] != "":
                    rej_seq_info[cds_id]["rejected_alternative_start"] = ",".join([
                        rej_seq_info[cds_id]["rejected_alternative_start"],
                        start])
                    rej_seq_info[cds_id]["rejected_alternative_end"] = ",".join([
                        rej_seq_info[cds_id]["rejected_alternative_end"],
                        end])
                else:
                    rej_seq_info[cds_id]["rejected_alternative_start"] = start
                    rej_seq_info[cds_id]["rejected_alternative_end"] = end
            else:
                rej_seq_info[cds_id]["start"] = start
                rej_seq_info[cds_id]["end"] = end
                rej_seq_info[cds_id]["strand"] = strand
    return conserved_seq, cons_seq_info, rej_seq_info


def check_pyl_proteins(
        potential_pyl_similarity_search, potential_pyl_seq,
        conserved_potential_pyl_sequences, conserved_potential_pyl_sequences_info,
        rejected_potential_pyl_sequences_info):
    """

    """
    cds_report = parse_similarity_search_report(
        potential_pyl_similarity_search)
    seq, cons_seq_info, rej_seq_info = extract_conserved_rejected_sequences(
        cds_report,
        potential_pyl_seq)

    with open(conserved_potential_pyl_sequences, "w") as conserved_seq_file:
        SeqIO.write(seq, conserved_seq_file, "fasta")

    cons_seq_info_df = pd.DataFrame(data=cons_seq_info).transpose()
    col = [
        "start",
        "end",
        "strand",
        "original_start",
        "original_end"]
    cons_seq_info_df = cons_seq_info_df[col]
    cons_seq_info_df.to_csv(conserved_potential_pyl_sequences_info)

    rej_seq_info_df = pd.DataFrame(data=rej_seq_info).transpose()
    col = [
        "start",
        "end",
        "strand",
        "rejected_alternative_start",
        "rejected_alternative_end",
        "comment"]
    rej_seq_info_df = rej_seq_info_df[col]
    rej_seq_info_df.to_csv(rejected_potential_pyl_sequences_info)
