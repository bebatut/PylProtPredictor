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
output_tag_ending_cds_filepath = os.path.join("tests","data","tag_ending_cds.csv")
output_cds_filepath = os.path.join("tests", "data", "cds.csv")
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

pred_cds, ordered_pred_cds, tag_ending_cds = predict_pyl_proteins.extract_predicted_cds(input_cds_filepath, "pred_cds", "tag_ending_cds",input_scaffold_genome_filepath)
#pot_pyl_cds = predict_pyl_proteins.extract_potential_pyl_cds(tag_ending_cds, pred_cds, input_scaffold_genome_filepath)


def test_get_string_seq():
    """Test get_string_seq"""
    s = predict_pyl_proteins.get_string_seq(genome[0])
    print(s)
    assert s.startswith("ATGCTCAGAAGAAAG")
    assert s.endswith("GATCGGGAAGGACGGAGC")


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


def test_extract_predicted_cds():
    """Test extract_predicted_cds"""
    assert filecmp.cmp("pred_cds", output_cds_filepath)
    assert "cds_1" in pred_cds
    assert "scaffold_0_3" in pred_cds
    assert "scaffold_0_1" in pred_cds
    assert "scaffold_117_243" in pred_cds
    assert "scaffold_1220_19" in pred_cds
    assert 'gi|477554117|gb|CP004049.1|' in ordered_pred_cds
    assert 'forward' in ordered_pred_cds['gi|477554117|gb|CP004049.1|']
    assert 'reverse' not in ordered_pred_cds['gi|477554117|gb|CP004049.1|']
    assert 'scaffold_0' in ordered_pred_cds
    assert 'forward' in ordered_pred_cds['scaffold_0']
    assert 'reverse' in ordered_pred_cds['scaffold_0']
    assert 'scaffold_117' in ordered_pred_cds
    assert 'forward' not in ordered_pred_cds['scaffold_117']
    assert 'reverse' in ordered_pred_cds['scaffold_117']
    assert 'scaffold_1220' in ordered_pred_cds
    assert 'forward' in ordered_pred_cds['scaffold_1220']
    assert 'reverse' not in ordered_pred_cds['scaffold_1220']
    assert filecmp.cmp("tag_ending_cds", output_tag_ending_cds_filepath)
    assert 'scaffold_117_243' in tag_ending_cds
    assert 'scaffold_1220_19' in tag_ending_cds
    os.remove('tag_ending_cds')
    os.remove('pred_cds')


def test_extract_seqs():
    """Test extract_seqs"""
    seq = predict_pyl_proteins.extract_seqs(input_cds_filepath)
    assert "gi|477554117|gb|CP004049.1|_1" in seq
    assert "scaffold_0_3" in seq
    assert len(seq["gi|477554117|gb|CP004049.1|_1"].seq) == 1245
    assert len(seq["scaffold_0_3"].seq) == 15


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


#def test_extract_potential_pyl_cds():
#    """Test extract_potential_pyl_cds"""
#    assert "scaffold_117_243" in pot_pyl_cds
#    assert len(pot_pyl_cds["scaffold_117_243"]["potential_seq"]) == 2
#    assert "scaffold_1220_19" in pot_pyl_cds
#    assert len(pot_pyl_cds["scaffold_1220_19"]["potential_seq"]) == 2


#def test_save_potential_pyl_cds():
#    """Test save_potential_pyl_cds"""
#    with open("log", 'w') as log_file:
#        predict_pyl_proteins.save_potential_pyl_cds(pot_pyl_cds, "pot_pyl_cds.fasta", log_file, "pot_pyl_cds")
#    assert filecmp.cmp("log", output_log_filepath)
#    assert filecmp.cmp("pot_pyl_cds.fasta", output_pot_pyl_cds_fasta_filepath)
#    assert filecmp.cmp("pot_pyl_cds", output_pot_pyl_cds_filepath)
#    os.remove('log')
#    os.remove('pot_pyl_cds.fasta')
#    os.remove('pot_pyl_cds')
