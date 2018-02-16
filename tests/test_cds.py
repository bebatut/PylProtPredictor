import pytest

from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq

from pylprotpredictor import cds


data_dir = Path("tests/data")
input_cds_filepath = str(data_dir / Path("cds.fasta"))
input_scaffold_genome_filepath = str(data_dir / Path("scaffold_genome.fasta"))
records = []
seqs = {}
cds_list = []
for record in SeqIO.parse(input_cds_filepath, "fasta"):
    records.append(record)
for record in SeqIO.parse(input_scaffold_genome_filepath, "fasta"):
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
    with pytest.raises(ValueError, match=r'Wrong strand_id: 2'):
        cds.transform_strand("2")


def test_find_stop_codon_pos_in_seq():
    """Test find_stop_codon_pos_in_seq"""
    seq = "OTGCAT*OTAAT*OO"
    c = cds.find_stop_codon_pos_in_seq(seq)
    assert len(c) == 2
    assert c[0] == 6
    assert c[1] == 12


def test_translate():
    """Test translate"""
    aa = cds.translate(records[1].seq)
    assert aa == "MOKO*"
    with pytest.raises(ValueError, match=r'Stop codon'):
        cds.translate(Seq("ATGTAGTGATGA"))


def test_init_from_record():
    """Test init_from_record function"""
    for record in SeqIO.parse(input_cds_filepath, "fasta"):
        cds_obj = cds.CDS()
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


def test_get_seqrecord():
    """Test get_seqrecord from CDS class"""
    seqrecord = cds_list[0].get_seqrecord()
    assert seqrecord.id == cds_list[0].get_id()
    assert seqrecord.seq.startswith("GTGGCGGATTCTATATTCAAA")


def test_get_translated_seq():
    """Test get_translated_seq from CDS class"""
    transl_seq = cds_list[0].get_translated_seq()
    assert transl_seq.id == cds_list[0].get_id()
    assert transl_seq.seq.startswith("VADSIFKRYLEKKNNL")


def test_set_id():
    """Test set_id from CDS class"""
    previous_id = cds_list[0].get_id()
    new_id = "test"
    cds_list[0].set_id(new_id)
    assert cds_list[0].get_id() == new_id
    cds_list[0].set_id(previous_id)
    assert cds_list[0].get_id() == previous_id


def test_set_origin_seq_id():
    """Test set_origin_seq_id from CDS class"""
    previous_id = cds_list[0].get_origin_seq_id()
    new_id = "test"
    cds_list[0].set_origin_seq_id(new_id)
    assert cds_list[0].get_origin_seq_id() == new_id
    cds_list[0].set_origin_seq_id(previous_id)
    assert cds_list[0].get_origin_seq_id() == previous_id


def test_set_start():
    """Test set_start from CDS class"""
    previous_start = cds_list[0].get_start()
    new_start = 10
    cds_list[0].set_start(new_start)
    assert cds_list[0].get_start() == new_start
    cds_list[0].set_start(previous_start)
    assert cds_list[0].get_start() == previous_start


def test_set_end():
    """Test set_end from CDS class"""
    previous_start = cds_list[0].get_end()
    new_start = 10
    cds_list[0].set_end(new_start)
    assert cds_list[0].get_end() == new_start
    cds_list[0].set_end(previous_start)
    assert cds_list[0].get_end() == previous_start


def test_set_strand():
    """Test set_strand from CDS class"""
    previous_strand = cds_list[0].get_strand()
    new_strand = "reverse"
    cds_list[0].set_strand(new_strand)
    assert cds_list[0].get_strand() == new_strand
    cds_list[0].set_strand(previous_strand)
    assert cds_list[0].get_strand() == previous_strand
    with pytest.raises(ValueError, match=r'Incorrect strand value: none'):
        cds_list[0].set_strand("none")


def test_set_seq():
    """Test set_seq from CDS class"""
    previous_seq = cds_list[0].get_seq()
    new_seq = "ATTGGCGGATGAC"
    cds_list[0].set_seq(new_seq)
    assert cds_list[0].get_seq() == new_seq
    cds_list[0].set_seq(previous_seq)
    assert cds_list[0].get_seq() == previous_seq


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


def test_get_origin_seq_size():
    """Test get_origin_seq_size from CDS class"""
    with pytest.raises(ValueError, match=r'No origin sequence provided'):
        cds_list[0].get_origin_seq_size()
    with pytest.raises(ValueError, match=r'No origin sequence provided'):
        cds_list[1].get_origin_seq_size() == 1
    with pytest.raises(ValueError, match=r'No origin sequence provided'):
        cds_list[2].get_origin_seq_size() == 1
    assert cds_list[3].get_origin_seq_size() == 245464
    assert cds_list[4].get_origin_seq_size() == 39692


def test_get_origin_seq_string():
    """Test get_origin_seq_string from CDS class"""
    with pytest.raises(ValueError, match=r'No origin sequence provided'):
        cds_list[0].get_origin_seq_string()
    assert cds_list[3].get_origin_seq_string().startswith("GTTTTGAGGT")
    assert cds_list[4].get_origin_seq_string().startswith("GATCAGAATT")


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


def test_has_next_cds_limit():
    """Test has_next_cds_limit from CDS class"""
    assert cds_list[0].has_next_cds_limit()
    assert not cds_list[1].has_next_cds_limit()


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
    with pytest.raises(ValueError, match=r'Incorrect order for scaffold_0_1'):
        cds_list[2].find_next_cds_limit(ordered_pred_cds, pred_cds)
    cds_list[2].set_order(2)
    with pytest.raises(ValueError, match=r'No origin sequence provided'):
        cds_list[2].find_next_cds_limit(ordered_pred_cds, pred_cds)
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
    """Test find_alternation_ends from CDS class"""
    with pytest.raises(ValueError, match=r'No origin sequence provided'):
        cds_list[0].find_alternative_ends()
    cds_list[1].set_origin_seq(seqs[cds_list[3].get_origin_seq_id()])
    cds_list[1].set_next_cds_limit(-1)
    with pytest.raises(ValueError, match=r'No next CDS limit provided'):
        cds_list[1].find_alternative_ends()
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


def test_add_alternative_cds():
    """Test get_alternative_cds from CDS class"""
    cds_list[0].add_alternative_cds(cds_list[1])
    assert len(cds_list[0].alternative_cds) == 1


def test_get_alternative_cds():
    """Test get_alternative_cds from CDS class"""
    alt_cds = cds_list[0].get_alternative_cds()
    assert len(alt_cds) == 1
    assert alt_cds[0].get_id() == cds_list[1].get_id()


def test_reset_alternative_cds():
    """Test reset_alternative_cds from CDS class"""
    cds_list[0].reset_alternative_cds()
    assert len(cds_list[0].get_alternative_cds()) == 0


def test_extract_possible_alternative_seq():
    """Test extract_possible_alternative_seq from CDS class"""
    with pytest.raises(ValueError, match=r'No origin sequence provided'):
        cds_list[0].extract_possible_alternative_seq()
    cds_list[3].extract_possible_alternative_seq()
    alt_cds = cds_list[3].get_alternative_cds()
    assert alt_cds[0].get_start() == 228199
    assert alt_cds[0].get_end() == 228441
    cds_list[4].extract_possible_alternative_seq()
    alt_cds = cds_list[4].get_alternative_cds()
    assert alt_cds[0].get_start() == 15451
    assert alt_cds[0].get_end() == 15957
    cds_list[2].set_alternative_ends([])
    cds_list[2].extract_possible_alternative_seq()
    assert len(cds_list[2].get_alternative_cds()) == 0


def test_get_alternative_start():
    """Test get_alternative_start from CDS class"""
    alt_start = cds_list[3].get_alternative_start()
    assert alt_start[0] == 228199
    alt_start = cds_list[4].get_alternative_start()
    assert alt_start[0] == 15451


def test_get_alternative_end():
    """Test get_alternative_end from CDS class"""
    alt_start = cds_list[3].get_alternative_end()
    assert alt_start[0] == 228441
    alt_start = cds_list[4].get_alternative_end()
    assert alt_start[0] == 15957


def test_get_translated_alternative_seq():
    """Test get_translated_alternative_seq from CDS class"""
    trans_alt_cds = cds_list[3].get_translated_alternative_seq()
    assert trans_alt_cds[0].seq.startswith("MKEIPFDDFINNISKAETP")
    print(trans_alt_cds[0].seq)
    assert trans_alt_cds[0].seq.count("O") == 1
    trans_alt_cds = cds_list[4].get_translated_alternative_seq()
    assert trans_alt_cds[0].seq.startswith("VANEEKKDFNAML")
    assert trans_alt_cds[0].seq.count("O") == 1


def test_export_description():
    """Test export_description from CDS class"""
    desc = cds_list[1].export_description()
    assert desc == "# origin_seq: scaffold_0 # strand: reverse # start: 1581 # end: 2864"


def test_is_potential_pyl_cds():
    """Test is_potential_pyl_cds from CDS class"""
    assert not cds_list[0].is_potential_pyl_cds()
    assert cds_list[3].is_potential_pyl_cds()
    assert cds_list[4].is_potential_pyl_cds()


def test_export_to_dict():
    """Test export_to_dict from CDS class"""
    d = cds_list[0].export_to_dict()
    cds_id = cds_list[0].get_id()
    assert cds_id in d
    assert d[cds_id]['id'] == cds_id
    assert d[cds_id]['origin_seq_id'] == cds_list[0].get_origin_seq_id()
    assert d[cds_id]['start'] == cds_list[0].get_start()
    assert d[cds_id]['end'] == cds_list[0].get_end()
    assert d[cds_id]['strand'] == cds_list[0].get_strand()
    assert d[cds_id]['seq'] == str(cds_list[0].get_seq())
    assert d[cds_id]['order'] == cds_list[0].get_order()
    assert d[cds_id]['next_cds_limit'] == cds_list[0].get_next_cds_limit()
    assert d[cds_id]['alternative_ends'] == cds_list[0].get_alternative_ends()
    assert d[cds_id]['alternative_cds'] == {}

    d = cds_list[4].export_to_dict()
    cds_id = cds_list[4].get_id()
    assert cds_id in d
    assert d[cds_id]['id'] == cds_id
    assert d[cds_id]['alternative_ends'] == cds_list[4].get_alternative_ends()
    assert 'scaffold_1220_19-1' in d[cds_id]['alternative_cds']
    assert d[cds_id]['alternative_cds']['scaffold_1220_19-1']['id'] == 'scaffold_1220_19-1'
    assert d[cds_id]['alternative_cds']['scaffold_1220_19-1']['alternative_ends'] == []


def test_init_from_dict():
    """Test init_from_dict from CDS class"""
    d = cds_list[0].export_to_dict()
    for k in d:
        new_cds = cds.CDS()
        new_cds.init_from_dict(d[k])
        assert new_cds.get_id() == cds_list[0].get_id()
        assert new_cds.get_origin_seq_id() == cds_list[0].get_origin_seq_id()
        assert new_cds.get_start() == cds_list[0].get_start()
        assert new_cds.get_end() == cds_list[0].get_end()
        assert new_cds.get_strand() == cds_list[0].get_strand()
        assert str(new_cds.get_seq()) == str(cds_list[0].get_seq())
        assert new_cds.get_order() == cds_list[0].get_order()
        assert new_cds.get_next_cds_limit() == cds_list[0].get_next_cds_limit()
        assert new_cds.has_alternative_ends()
        assert new_cds.get_alternative_ends() == cds_list[0].get_alternative_ends()
        assert not new_cds.is_potential_pyl_cds()

    d = cds_list[4].export_to_dict()
    for k in d:
        new_cds = cds.CDS()
        new_cds.init_from_dict(d[k])
        assert new_cds.get_id() == cds_list[4].get_id()
        assert new_cds.has_alternative_ends()
        assert new_cds.get_alternative_ends() == cds_list[4].get_alternative_ends()
        assert new_cds.is_potential_pyl_cds()


def test_set_evalue():
    """Test set_evalue from CDS class"""
    cds_list[0].set_evalue(0.1)
    assert cds_list[0].evalue == 0.1


def test_get_evalue():
    """Test get_evalue from CDS class"""
    assert cds_list[0].get_evalue() == 0.1
    assert cds_list[2].get_evalue() == 10


def test_set_conserved_cds():
    """Test set_conserved_cds from CDS class"""
    cds_list[0].set_conserved_cds(cds_list[1])
    assert cds_list[0].conserved_cds == cds_list[1]


def test_get_conserved_cds():
    """Test get_conserved_cds from CDS class"""
    assert cds_list[0].get_conserved_cds() == cds_list[1]
    assert cds_list[2].get_conserved_cds() is None


def test_add_rejected_cds():
    """Test get_rejected_cds from CDS class"""
    cds_list[0].add_rejected_cds(cds_list[1])
    assert len(cds_list[0].rejected_cds) == 1


def test_get_rejected_cds():
    """Test get_rejected_cds from CDS class"""
    alt_cds = cds_list[0].get_rejected_cds()
    assert len(alt_cds) == 1
    assert alt_cds[0].get_id() == cds_list[1].get_id()


def test_reset_rejected_cds():
    """Test reset_alternative_cds from CDS class"""
    cds_list[0].reset_rejected_cds()
    assert len(cds_list[0].get_rejected_cds()) == 0


def test_add_evalue():
    """Test add_evalue from CDS class"""
    cds_list[0].add_evalue(cds_list[0].get_id(), 0.02)
    assert cds_list[0].get_evalue() == 0.02
    cds_list[4].add_evalue("scaffold_1220_19-1", 0.02)
    for alt_cds in cds_list[4].get_alternative_cds():
        assert alt_cds.get_evalue() == 0.02


def test_identify_cons_rej_cds():
    """Test identify_cons_rej_cds from CDS class"""
    cds_list[3].identify_cons_rej_cds()
    assert cds_list[3].get_conserved_cds().get_id() == cds_list[3].get_id()
    assert len(cds_list[3].get_rejected_cds()) == 1
    cds_list[4].identify_cons_rej_cds()
    assert cds_list[4].get_conserved_cds().get_id() == "scaffold_1220_19-1"
    assert cds_list[4].get_rejected_cds()[0].get_id() == cds_list[4].get_id()
