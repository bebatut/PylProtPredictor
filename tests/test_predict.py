import filecmp
import os

from pathlib import Path

from Bio import SeqIO

from pylprotpredictor import predict


data_dir = Path("tests/data")
input_cds_filepath = data_dir / Path("cds.fasta")
input_genome_filepath = data_dir / Path("genome.fasta")
input_scaffold_genome_filepath = data_dir / Path("scaffold_genome.fasta")
output_tag_ending_cds_filepath = data_dir / Path("tag_ending_cds.csv")
output_cds_filepath = data_dir / Path("cds.csv")
output_log_filepath = data_dir / Path("log")
output_pot_pyl_cds_filepath = data_dir / Path("pot_pyl_cds")
output_pot_pyl_cds_fasta_filepath = data_dir / Path("pot_pyl_cds.fasta")
output_pred_cds_obj_filepath = data_dir / Path("pred_cds_obj.json")
descs = []
seqs = []
for record in SeqIO.parse(input_cds_filepath.as_posix(), "fasta"):
    descs.append(record.description)
    seqs.append(record.seq)
genome = []
for record in SeqIO.parse(input_genome_filepath.as_posix(), "fasta"):
    genome.append(record)

pred_cds = predict.extract_predicted_cds(
    input_cds_filepath,
    Path("pred_cds"),
    Path("tag_ending_cds"),
    input_scaffold_genome_filepath)


def test_extract_predicted_cds():
    """Test extract_predicted_cds"""
    assert filecmp.cmp("pred_cds", output_cds_filepath.as_posix())
    assert "cds_1" in pred_cds
    assert "scaffold_0_3" in pred_cds
    assert "scaffold_0_1" in pred_cds
    assert "scaffold_117_243" in pred_cds
    assert "scaffold_1220_19" in pred_cds
    assert filecmp.cmp("tag_ending_cds", output_tag_ending_cds_filepath.as_posix())
    os.remove('tag_ending_cds')
    os.remove('pred_cds')


def test_extract_seqs():
    """Test extract_seqs"""
    seq = predict.extract_seqs(input_cds_filepath)
    assert "gi|477554117|gb|CP004049.1|_1" in seq
    assert "scaffold_0_3" in seq
    assert len(seq["gi|477554117|gb|CP004049.1|_1"].seq) == 1245
    assert len(seq["scaffold_0_3"].seq) == 15


def test_extract_potential_pyl_cds():
    """Test extract_potential_pyl_cds"""
    predict.extract_potential_pyl_cds(
        pred_cds,
        Path("pot_pyl_cds.fasta"),
        Path("pot_pyl_cds"),
        Path("pred_cds_obj"))
    assert os.stat("pot_pyl_cds.fasta").st_size == os.stat(output_pot_pyl_cds_fasta_filepath.as_posix()).st_size
    assert filecmp.cmp("pot_pyl_cds", output_pot_pyl_cds_filepath.as_posix())
    assert os.stat("pred_cds_obj").st_size == os.stat(output_pred_cds_obj_filepath.as_posix()).st_size
    os.remove('pot_pyl_cds.fasta')
    os.remove('pot_pyl_cds')
    os.remove('pred_cds_obj')


def test_predict_pyl_proteins():
    """Test predict_pyl_proteins"""
    predict.predict_pyl_proteins(
        input_scaffold_genome_filepath,
        input_cds_filepath,
        Path("pot_pyl_cds.fasta"),
        Path("log"),
        Path("pred_cds"),
        Path("tag_ending_cds"),
        Path("pot_pyl_cds"),
        Path("pred_cds_obj"))
    assert filecmp.cmp("pred_cds", output_cds_filepath.as_posix())
    assert filecmp.cmp("tag_ending_cds", output_tag_ending_cds_filepath.as_posix())
    assert os.stat("pot_pyl_cds.fasta").st_size == os.stat(output_pot_pyl_cds_fasta_filepath.as_posix()).st_size
    assert filecmp.cmp("pot_pyl_cds", output_pot_pyl_cds_filepath.as_posix())
    assert os.stat("pred_cds_obj").st_size == os.stat(output_pred_cds_obj_filepath.as_posix()).st_size
    os.remove('pred_cds')
    os.remove('tag_ending_cds')
    os.remove('pot_pyl_cds.fasta')
    os.remove('pot_pyl_cds')
    # os.remove('pred_cds_obj')
