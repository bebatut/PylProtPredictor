#!/usr/bin/env python

import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from src import cds


input_cds_filepath = os.path.join("tests","data","cds.fasta")
input_scaffold_genome_filepath = os.path.join("tests","data","scaffold_genome.fasta")
records = []
cds_list = []
for record in SeqIO.parse(input_cds_filepath,"fasta"):
    records.append(record)seqs = {}
for record in SeqIO.parse(input_scaffold_genome_filepath,"fasta"):
    seqs[record.id] = record


def test_extract_seq_desc():
    """Test extract_seq_desc function"""
    # first test
    seq_id, origin_seq, start, end, strand = cds.extract_seq_desc(records[0].description)
    assert seq_id == "cds_1"
    assert origin_seq == "gi|477554117|gb|CP004049.1|"
    assert start == 694
    assert end == 1938
    assert strand == "forward"
    # second test (reverse)
    seq_id, origin_seq, start, end, strand = cds.extract_seq_desc(records[1].description)
    assert seq_id == "scaffold_0_3"
    assert origin_seq == "scaffold_0"
    assert start == 1581
    assert end == 2864
    assert strand == "reverse"


def test_transform_strand():
    """Test transform_strand function"""
    res = cds.transform_strand("-1")
    assert res == "reverse"
    print("Wrong strand_id: {}".format(res))
    res = cds.transform_strand("1")
    assert res == "forward"


def test_init_from_record():
    """Test init_from_record function"""
    for record in SeqIO.parse(input_cds_filepath,"fasta"):
        cds_obj = cds()
        cds_obj.init_from_record(record)
        cds_list.append(cds_obj)
    assert cds_list[0].id == "cds_1"


def test_get_id():
    """Test get_id from CDS class"""
    assert cds_list[0].get_id() == "cds_1"
    assert cds_list[1].get_id() == "scaffold_0_3"
    assert cds_list[2].get_id() == "scaffold_0_1"
    assert cds_list[3].get_id() == "scaffold_117_243"
    assert cds_list[4].get_id() == "scaffold_1220_19"


def test_get_origin_seq_id():
    """Test get_origin from CDS class"""
    assert cds_list[0].get_origin_seq_id() == "gi|477554117|gb|CP004049.1|"
    assert cds_list[1].get_origin_seq_id() == "scaffold_0"
    assert cds_list[2].get_origin_seq_id() == "scaffold_0"
    assert cds_list[3].get_origin_seq_id() == "scaffold_117"
    assert cds_list[4].get_origin_seq_id() == "scaffold_1220"


def test_get_start():
    """Test get_start from CDS class"""
    assert cds_list[0].get_start() == 694
    assert cds_list[1].get_start() == 1581
    assert cds_list[2].get_start() == 1
    assert cds_list[3].get_start() == 228271
    assert cds_list[4].get_start() == 15451


def test_get_end():
    """Test get_end from CDS class"""
    assert cds_list[0].get_end() == 1938
    assert cds_list[1].get_end() == 2864
    assert cds_list[2].get_end() == 543
    assert cds_list[3].get_end() == 228441
    assert cds_list[4].get_end() == 15930


def test_get_strand():
    """Test get_strand from CDS class"""
    assert cds_list[0].get_strand() == "forward"
    assert cds_list[1].get_strand() == "reverse"
    assert cds_list[2].get_strand() == "forward"
    assert cds_list[3].get_strand() == "reverse"
    assert cds_list[4].get_strand() == "forward"


def test_get_seq():
    """Test get_seq from CDS class"""
    assert cds_list[0].get_seq().startswith("GTGGCGGATTCTATATTCAAA")
    assert cds_list[1].get_seq() == "ATGTAGAAATAGTAA"
    assert cds_list[2].get_seq().startswith("GATTCAATTTCAATAGCA")
    assert cds_list[3].get_seq().startswith("ATGAAAGAAATACC")
    assert cds_list[4].get_seq().startswith("GTGGCAAATGAA")


def test_is_tag_ending_seq():
    """Test is_tag_ending from CDS class"""
    assert not cds_list[0].is_tag_ending_seq()
    assert not cds_list[1].is_tag_ending_seq()
    assert not cds_list[2].is_tag_ending_seq()
    assert cds_list[3].is_tag_ending_seq()
    assert cds_list[4].is_tag_ending_seq()


def test_set_order():
    """Test set_order from CDS class"""
    cds_list[0].set_order(2)
    assert cds_list[0].order == 2


def test_get_order():
    """Test get_order from CDS class"""
    assert cds_list[0].get_order() == 2


def test_set_origin_seq():
    """Test set_origin_seq from CDS class"""
    origin_seq = cds_list[3].get_origin_seq_id()
    cds_list[3].set_origin_seq(seqs[origin_seq])
    assert cds_list[3].origin_seq is not None
    origin_seq = cds_list[4].get_origin_seq_id()
    cds_list[4].set_origin_seq(seqs[origin_seq])
    assert cds_list[4].origin_seq is not None
    assert cds_list[0].origin_seq is None


def test_get_origin_seq():
    """Test get_origin_seq from CDS class"""
    origin_seq = cds_list[3].get_origin_seq()
    assert origin_seq.description == "scaffold_117"
    assert origin_seq.seq.startswith("GTTTTGAGGT")
    origin_seq = cds_list[4].get_origin_seq()
    assert origin_seq.description == "scaffold_1220"
    assert origin_seq.seq.startswith("GATCAGAATT")
    origin_seq = cds_list[0].get_origin_seq()
    assert origin_seq is None


def test_is_reverse_strand():
    """Test is_reverse_strand from CDS class"""
    assert not cds_list[0].is_reverse_strand()
    assert cds_list[1].is_reverse_strand()
    assert not cds_list[2].is_reverse_strand()
    assert cds_list[3].is_reverse_strand()
    assert not cds_list[4].is_reverse_strand()


def test_set_next_cds_limit():
    """Test set_next_cds_limit from CDS class"""
    cds_list[0].set_next_cds_limit(2)
    assert cds_list[0].next_cds_limit == 2


def test_get_next_cds_limit():
    """Test get_next_cds_limit from CDS class"""
    assert cds_list[0].get_next_cds_limit() == 2


def test_find_next_cds_limit():
    """Test find_next_cds_limit from CDS class"""
    ordered_pred_cds = ['cds_1', 'scaffold_0_3', 'scaffold_0_1', 'scaffold_117_243', 'scaffold_1220_19']
    pred_cds = {
        'cds_1': cds_list[0],
        'scaffold_0_3': cds_list[1],
        'scaffold_0_1': cds_list[2],
        'scaffold_117_243': cds_list[3],
        'scaffold_1220_19': cds_list[4]
    }
    cds_list[2].set_order(2)
    cds_list[2].set_origin_seq(seqs[cds_list[3].get_origin_seq_id()])
    cds_list[2].find_next_cds_limit(ordered_pred_cds, pred_cds)
    assert cds_list[2].get_next_cds_limit() == 228441
    cds_list[3].set_order(3)
    cds_list[3].find_next_cds_limit(ordered_pred_cds, pred_cds)
    assert cds_list[3].get_next_cds_limit() == 245464
    cds_list[4].set_order(4)
    cds_list[4].find_next_cds_limit(ordered_pred_cds, pred_cds)
    assert cds_list[4].get_next_cds_limit() == len(cds_list[4].get_origin_seq().seq)


def test_set_alternative_ends():
    """Test set_alternative_ends from CDS class"""
    cds_list[0].set_alternative_ends([2, 25533])
    assert len(cds_list[0].alternative_ends) == 2


def test_get_alternative_ends():
    """Test get_alternative_ends from CDS class"""
    alt_end = cds_list[0].get_alternative_ends()
    assert len(alt_end) == 2
    assert alt_end[0] == 2
    assert alt_end[1] == 25533


def test_has_alternative_ends():
    """Test has_alternative_ends from CDS class"""
    assert cds_list[0].has_alternative_ends()
    assert not cds_list[1].has_alternative_ends()


def test_find_alternative_ends():
    """Test find_alternation_ends"""
    cds_list[2].find_alternative_ends()
    alt_end = cds_list[2].get_alternative_ends()
    assert len(alt_end) == 1
    assert alt_end[0] == 552
    cds_list[3].find_alternative_ends()
    alt_end = cds_list[3].get_alternative_ends()
    assert len(alt_end) == 1
    assert alt_end[0] == 17266
    cds_list[4].find_alternative_ends()
    alt_end = cds_list[4].get_alternative_ends()
    assert len(alt_end) == 1
    assert alt_end[0] == 15957


