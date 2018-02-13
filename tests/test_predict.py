import filecmp
import os

from pathlib import Path

from Bio import SeqIO

from pylprotpredictor import predict


data_dir = Path("tests/data")
input_cds_filepath = str(data_dir / Path("cds.fasta"))
input_genome_filepath = str(data_dir / Path("genome.fasta"))
input_scaffold_genome_filepath = str(data_dir / Path("scaffold_genome.fasta"))
output_tag_ending_cds_filepath = str(data_dir / Path("tag_ending_cds.csv"))
output_cds_filepath = str(data_dir / Path("cds.csv"))
output_log_filepath = str(data_dir / Path("log"))
output_pot_pyl_cds_filepath = str(data_dir / Path("pot_pyl_cds"))
output_pot_pyl_cds_fasta_filepath = str(data_dir / Path("pot_pyl_cds.fasta"))
output_pot_pyl_obj_filepath = str(data_dir / Path("pot_pyl_obj.json"))
descs = []
seqs = []
for record in SeqIO.parse(input_cds_filepath, "fasta"):
    descs.append(record.description)
    seqs.append(record.seq)
genome = []
for record in SeqIO.parse(input_genome_filepath, "fasta"):
    genome.append(record)

pred_cds, ordered_pred_cds, tag_ending_cds = predict.extract_predicted_cds(
    input_cds_filepath,
    "pred_cds",
    "tag_ending_cds",
    input_scaffold_genome_filepath)


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
    seq = predict.extract_seqs(input_cds_filepath)
    assert "gi|477554117|gb|CP004049.1|_1" in seq
    assert "scaffold_0_3" in seq
    assert len(seq["gi|477554117|gb|CP004049.1|_1"].seq) == 1245
    assert len(seq["scaffold_0_3"].seq) == 15


def test_extract_potential_pyl_cds():
    """Test extract_potential_pyl_cds"""
    predict.extract_potential_pyl_cds(
        tag_ending_cds,
        pred_cds,
        ordered_pred_cds,
        "pot_pyl_cds.fasta",
        "pot_pyl_cds",
        "pot_pyl_obj")
    assert filecmp.cmp("pot_pyl_cds.fasta", output_pot_pyl_cds_fasta_filepath)
    assert filecmp.cmp("pot_pyl_cds", output_pot_pyl_cds_filepath)
    assert os.stat("pot_pyl_obj").st_size == os.stat(output_pot_pyl_obj_filepath).st_size
    os.remove('pot_pyl_cds.fasta')
    os.remove('pot_pyl_cds')
    os.remove('pot_pyl_obj')


def test_predict_pyl_proteins():
    """Test predict_pyl_proteins"""
    predict.predict_pyl_proteins(
        input_scaffold_genome_filepath,
        input_cds_filepath,
        "pot_pyl_cds.fasta",
        "log",
        "pred_cds",
        "tag_ending_cds",
        "pot_pyl_cds",
        "pot_pyl_obj")
    assert filecmp.cmp("pred_cds", output_cds_filepath)
    assert filecmp.cmp("tag_ending_cds", output_tag_ending_cds_filepath)
    assert filecmp.cmp("pot_pyl_cds.fasta", output_pot_pyl_cds_fasta_filepath)
    assert filecmp.cmp("pot_pyl_cds", output_pot_pyl_cds_filepath)
    assert os.stat("pot_pyl_obj").st_size == os.stat(output_pot_pyl_obj_filepath).st_size
    os.remove('pred_cds')
    os.remove('tag_ending_cds')
    os.remove('pot_pyl_cds.fasta')
    os.remove('pot_pyl_cds')
    os.remove('pot_pyl_obj')
