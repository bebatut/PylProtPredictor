import os

from pathlib import Path

from pylprotpredictor import write_report


data_dir = Path("tests/data")
pred_cds_filepath = str(data_dir / Path("cds.csv"))
tag_ending_cds_filepath = str(data_dir / Path("tag_ending_cds.csv"))
pot_pyl_protein_filepath = str(data_dir / Path("pot_pyl_cds"))
final_cds_filepath = str(data_dir / Path("final_info.csv"))
report_filepath = str(data_dir / Path("report.html"))


def test_extract_row_number():
    """Test extract_row_number"""
    assert write_report.extract_row_number(pred_cds_filepath) == 5
    assert write_report.extract_row_number(tag_ending_cds_filepath) == 2
    assert write_report.extract_row_number(pot_pyl_protein_filepath) == 2
    assert write_report.extract_row_number(final_cds_filepath) == 2


def test_write_report():
    """Test write_report"""
    write_report.write_report(
        pred_cds_filepath,
        tag_ending_cds_filepath,
        pot_pyl_protein_filepath,
        final_cds_filepath,
        "report.html")
    assert round(os.stat("report.html").st_size, -1) == round(os.stat(report_filepath).st_size, -1)
    os.remove("report.html")
