import filecmp
import os
import pytest

from pathlib import Path

from pylprotpredictor import check


data_dir = Path("tests/data")
pot_pyl_cds_filepath = str(data_dir / Path("pot_pyl_obj.json"))
sim_search_filepath = str(data_dir / Path("similarity_search.txt"))
cons_pot_pyl_seq_filepath = str(data_dir / Path("cons_pot_pyl_seq.fasta"))
info_filepath = str(data_dir / Path("final_info.csv"))

pot_pyl_cds = check.import_cds(pot_pyl_cds_filepath)


def test_import_cds():
    """Test import_cds function"""
    assert "scaffold_117_243" in pot_pyl_cds
    assert pot_pyl_cds["scaffold_117_243"].get_id() == "scaffold_117_243"


def test_get_cds_obj():
    """Test get_cds_obj function"""
    cds_obj = check.get_cds_obj("scaffold_117_243", pot_pyl_cds)
    assert cds_obj.get_id() == "scaffold_117_243"
    with pytest.raises(ValueError, match=r'CDS not found for scaffold_117'):
        check.get_cds_obj("scaffold_117", pot_pyl_cds)


def test_parse_similarity_search_report():
    """Test parse_similarity_search_report function"""
    check.parse_similarity_search_report(sim_search_filepath, pot_pyl_cds)
    assert pot_pyl_cds["scaffold_117_243"].get_alignments()[0].get_evalue() == 7.000000000000001e-54
    assert pot_pyl_cds["scaffold_117_243"].get_alternative_cds()[0].get_alignments()[0].get_evalue() == 3.999999999999999e-131


def test_extract_correct_cds():
    """Test extract_correct_cds function"""
    check.extract_correct_cds(pot_pyl_cds, "cons_pot_pyl_seq", "info")
    assert filecmp.cmp("cons_pot_pyl_seq", cons_pot_pyl_seq_filepath)
    assert filecmp.cmp("info", info_filepath)
    os.remove("cons_pot_pyl_seq")
    os.remove("info")


def test_check_pyl_proteins():
    """Test check_pyl_proteins function"""
    check.check_pyl_proteins(sim_search_filepath, pot_pyl_cds_filepath, "cons_pot_pyl_seq", "info")
    assert filecmp.cmp("cons_pot_pyl_seq", cons_pot_pyl_seq_filepath)
    assert filecmp.cmp("info", info_filepath)
    os.remove("cons_pot_pyl_seq")
    os.remove("info")
