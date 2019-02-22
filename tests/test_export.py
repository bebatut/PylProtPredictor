import filecmp
import os

from pathlib import Path

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from pylprotpredictor import export


data_dir = Path("tests/data")
pot_pyl_cds = data_dir / Path("pot_pyl_cds.fasta")
cds = data_dir / Path("cds.csv")
light_cds = data_dir / Path("light_cds.csv")
test = data_dir / Path("test.json")

info = {
    'cds_1': {'start': 694, 'end': 1938, 'strand': 'forward', 'origin_seq': 'gi|477554117|gb|CP004049.1|', 'stop_codon': 'TGA'},
    'scaffold_0_3': {'start': 1581, 'end': 2864, 'strand': 'reverse', 'origin_seq': 'scaffold_0', 'stop_codon': 'TAA'},
    'scaffold_0_1': {'start': 1, 'end': 543, 'strand': 'forward', 'origin_seq': 'scaffold_0', 'stop_codon': 'TGA'},
    'scaffold_117_243': {'start': 228271, 'end': 228441, 'strand': 'reverse', 'origin_seq': 'scaffold_117', 'stop_codon': 'TAG'},
    'scaffold_1220_19': {'start': 15451, 'end': 15930, 'strand': 'forward', 'origin_seq': 'scaffold_1220', 'stop_codon': 'TAG'}
}


def test_export_csv():
    """Test export_csv function"""
    # Export all columns
    tmp = Path("tmp")
    export.export_csv(info, tmp, ["start", "end", "strand", "origin_seq", 'stop_codon'])
    assert filecmp.cmp(tmp.as_posix(), cds.as_posix())
    os.remove(tmp.as_posix())
    # Export only two columns
    export.export_csv(info, tmp, ["start", "end"])
    assert filecmp.cmp(tmp.as_posix(), light_cds.as_posix())
    os.remove(tmp.as_posix())


def test_export_fasta():
    """Test export_fasta function"""
    tmp = Path("tmp")
    seq = [
        SeqRecord(
            Seq("VANEEKKDFNAMLRKNTDMPKTQIVTDESIIKRYGGERMFFAPPLAYDELMKKVPHGKVVTAEKIREYLAEKNCADFTDPMTAGLFISIAAWASHQREEDITPYWRTLKTDGELNAKYPGGIEAQKKMLEEEGHVIIQKGRKNIRFFVKDYENVLFDLH"),
            id="scaffold_1220_19",
            description="# origin_seq: scaffold_1220 # strand: forward # start: 15451 # end: 15930"),
        SeqRecord(
            Seq("VANEEKKDFNAMLRKNTDMPKTQIVTDESIIKRYGGERMFFAPPLAYDELMKKVPHGKVVTAEKIREYLAEKNCADFTDPMTAGLFISIAAWASHQREEDITPYWRTLKTDGELNAKYPGGIEAQKKMLEEEGHVIIQKGRKNIRFFVKDYENVLFDLHOSDLVKLFL"),
            id="scaffold_1220_19-1",
            description="# origin_seq: scaffold_1220 # strand: forward # start: 15451 # end: 15957"),
        SeqRecord(
            Seq("MKEIPFDDFINNISKAETPFAVKIDNMLLQGNCLREIKEGASGYSVRMKRHESISW"),
            id="scaffold_117_243",
            description="# origin_seq: scaffold_117 # strand: reverse # start: 228271 # end: 228441"),
        SeqRecord(
            Seq("MKEIPFDDFINNISKAETPFAVKIDNMLLQGNCLREIKEGASGYSVRMKRHESISWOLPDYKSNKGVFKSELISVLIFQV"),
            id="scaffold_117_243-1",
            description="# origin_seq: scaffold_117 # strand: reverse # start: 228199 # end: 228441")]
    export.export_fasta(seq, tmp)
    assert filecmp.cmp(tmp.as_posix(), pot_pyl_cds.as_posix())
    os.remove(tmp.as_posix())


def test_export_json():
    """Test export_json function"""
    tmp = Path("tmp")
    export.export_json(info, tmp)
    assert os.stat(tmp.as_posix()).st_size == os.stat(test.as_posix()).st_size
    os.remove(tmp.as_posix())


def test_import_json():
    """Test export_json function"""
    d = export.import_json(test)
    for s in info:
        assert s in d
        for t in info[s]:
            assert t in d[s]
