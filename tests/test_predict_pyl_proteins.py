#!/usr/bin/env python

from src import predict_pyl_proteins
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import filecmp


input_cds_filepath = os.path.join("tests","data","cds.fasta")
input_genome_filepath = os.path.join("tests","data","genome.fasta")
input_scaffold_genome_filepath = os.path.join("tests","data","scaffold_genome.fasta")
output_cds_filepath = os.path.join("tests","data","cds.csv")
output_tag_ending_cds_filepath = os.path.join("tests","data","tag_ending_cds.csv")
output_log_filepath = os.path.join("tests","data","log")
output_pot_pyl_cds_filepath = os.path.join("tests","data","pot_pyl_cds")
output_pot_pyl_cds_fasta_filepath = os.path.join("tests","data","pot_pyl_cds.fasta")
descs = []
seqs = []
for record in SeqIO.parse(input_cds_filepath,"fasta"):
    descs.append(record.description)
    seqs.append(record.seq)
genome = []
for record in SeqIO.parse(input_genome_filepath,"fasta"):
    genome.append(record)

pred_cds, cds_nb = predict_pyl_proteins.extract_predicted_cds(input_cds_filepath, "pred_cds")
tag_ending_cds, tag_ending_cds_nb = predict_pyl_proteins.identify_tag_ending_cds(pred_cds, "tag_ending_cds")
pot_pyl_cds = predict_pyl_proteins.extract_potential_pyl_cds(tag_ending_cds, pred_cds, input_scaffold_genome_filepath)


def test_transform_strand():
    """Test transform_strand function"""
    res = predict_pyl_proteins.transform_strand("-1")
    assert res == "reverse"
    print("Wrong strand_id: {}".format(res))
    res = predict_pyl_proteins.transform_strand("1")
    assert res == "forward"


def test_get_string_seq():
    """Test get_string_seq"""
    s = predict_pyl_proteins.get_string_seq(genome[0])
    print(s)
    assert s.startswith("ATGCTCAGAAGAAAG")
    assert s.endswith("GATCGGGAAGGACGGAGC")


def test_extract_seq_desc():
    """Test extract_seq_desc function"""
    # first test
    seq_id, origin_seq, start, end, strand = predict_pyl_proteins.extract_seq_desc(descs[0])
    assert seq_id == "cds_1"
    assert origin_seq == "gi|477554117|gb|CP004049.1|"
    assert start == 694
    assert end == 1938
    assert strand == "forward"
    # second test (reverse)
    seq_id, origin_seq, start, end, strand = predict_pyl_proteins.extract_seq_desc(descs[1])
    assert seq_id == "scaffold_0_3"
    assert origin_seq == "scaffold_0"
    assert start == 1581
    assert end == 2864
    assert strand == "reverse"


def test_translate():
    """Test translate"""
    aa = predict_pyl_proteins.translate(seqs[1])
    assert aa == "MOKO*"


def test_find_stop_codon_pos_in_seq():
    """Test find_stop_codon_pos_in_seq"""
    seq = "OTGCAT*OTAAT*OO"
    c = predict_pyl_proteins.find_stop_codon_pos_in_seq(seq)
    assert len(c) == 2
    assert c[0] == 6
    assert c[1] == 12


def test_export_csv():
    """Test export_csv function"""
    pred_cds_info = {}
    for desc in descs:
        seq_id, origin_seq, start, end, strand = predict_pyl_proteins.extract_seq_desc(desc)
        pred_cds_info[seq_id] = {"start": start, "end": end, "strand": strand, "origin_seq": origin_seq}
    predict_pyl_proteins.export_csv(pred_cds_info, "tmp", ["start", "end", "strand", "origin_seq"])
    assert filecmp.cmp("tmp", output_cds_filepath)
    os.remove('tmp')


def test_extract_predicted_cds():
    """Test extract_predicted_cds"""
    assert filecmp.cmp("pred_cds", output_cds_filepath)
    assert 'gi|477554117|gb|CP004049.1|' in pred_cds
    assert 'forward' in pred_cds['gi|477554117|gb|CP004049.1|']
    assert 'reverse' not in pred_cds['gi|477554117|gb|CP004049.1|']
    assert 'scaffold_0' in pred_cds
    assert 'forward' in pred_cds['scaffold_0']
    assert 'reverse' in pred_cds['scaffold_0']
    assert 'scaffold_117' in pred_cds
    assert 'forward' not in pred_cds['scaffold_117']
    assert 'reverse' in pred_cds['scaffold_117']
    assert 'scaffold_1220' in pred_cds
    assert 'forward' in pred_cds['scaffold_1220']
    assert 'reverse' not in pred_cds['scaffold_1220']
    assert cds_nb == 5
    os.remove('pred_cds')


def test_extract_seqs():
    """Test extract_seqs"""
    seq = predict_pyl_proteins.extract_seqs(input_cds_filepath)
    assert "gi|477554117|gb|CP004049.1|_1" in seq
    assert "scaffold_0_3" in seq
    assert len(seq["gi|477554117|gb|CP004049.1|_1"]["rev_comp_genome"]) == 1245
    assert len(seq["scaffold_0_3"]["genome"]) == 15


def test_find_alternative_ends():
    """Test find_alternation_ends"""
    ends = predict_pyl_proteins.find_alternative_ends(300, predict_pyl_proteins.get_string_seq(genome[0]), 500)
    assert ends == [324, 354]
    ends = predict_pyl_proteins.find_alternative_ends(300, predict_pyl_proteins.get_string_seq(genome[0]), 350)
    assert ends == [324]


def test_identify_tag_ending_cds():
    """Test identify_tag_ending_proteins"""
    assert filecmp.cmp("tag_ending_cds", output_tag_ending_cds_filepath)
    assert 'scaffold_117' in tag_ending_cds
    assert 'forward' not in tag_ending_cds['scaffold_117']
    assert 'reverse' in tag_ending_cds['scaffold_117']
    assert 'scaffold_117_243' in tag_ending_cds['scaffold_117']['reverse']
    assert 'scaffold_1220' in tag_ending_cds
    assert 'forward' in tag_ending_cds['scaffold_1220']
    assert 'reverse' not in tag_ending_cds['scaffold_1220']
    assert 'scaffold_1220_19' in tag_ending_cds['scaffold_1220']['forward']
    assert tag_ending_cds_nb == 2
    os.remove('tag_ending_cds')


def test_extract_possible_seq():
    """Test extract_possible_seq"""
    pot_pyl_prot_def = predict_pyl_proteins.extract_possible_seq(
        200,
        300,
        predict_pyl_proteins.find_alternative_ends(300, genome[0].seq, 500),
        "forward",
        genome[0].id,
        predict_pyl_proteins.get_string_seq(genome[0]))
    assert pot_pyl_prot_def['strand'] == 'forward'
    assert pot_pyl_prot_def['origin_seq'] == 'test'
    assert len(pot_pyl_prot_def['potential_seq']) == 3


def test_extract_potential_pyl_cds():
    """Test extract_potential_pyl_cds"""
    assert "scaffold_117_243" in pot_pyl_cds
    assert len(pot_pyl_cds["scaffold_117_243"]["potential_seq"]) == 2
    assert "scaffold_1220_19" in pot_pyl_cds
    assert len(pot_pyl_cds["scaffold_1220_19"]["potential_seq"]) == 2


def test_save_potential_pyl_cds():
    """Test save_potential_pyl_cds"""
    with open("log", 'w') as log_file:
        predict_pyl_proteins.save_potential_pyl_cds(pot_pyl_cds, "pot_pyl_cds.fasta", log_file, "pot_pyl_cds")
    assert filecmp.cmp("log", output_log_filepath)
    assert filecmp.cmp("pot_pyl_cds.fasta", output_pot_pyl_cds_fasta_filepath)
    assert filecmp.cmp("pot_pyl_cds", output_pot_pyl_cds_filepath)
    os.remove('log')
    os.remove('pot_pyl_cds.fasta')
    os.remove('pot_pyl_cds')
