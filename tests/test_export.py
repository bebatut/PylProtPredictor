#!/usr/bin/env python

from src import export
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import filecmp
import os


def test_export_csv():
    """Test export_csv function"""
    info = {
        'cds_1': {'start': 694, 'end': 1938, 'strand': 'forward', 'origin_seq': 'gi|477554117|gb|CP004049.1|'},
        'scaffold_0_1': {'start': 1, 'end': 543, 'strand': 'forward', 'origin_seq': 'scaffold_0'},
        'scaffold_0_3': {'start': 1581, 'end': 2864, 'strand': 'reverse', 'origin_seq': 'scaffold_0'},
        'scaffold_117_243': {'start': 228271, 'end': 228441, 'strand': 'reverse', 'origin_seq': 'scaffold_117'},
        'scaffold_1220_19': {'start': 15451, 'end': 15930, 'strand': 'forward', 'origin_seq': 'scaffold_1220'}
    }
    # Export all columns
    export.export_csv(info, "tmp", ["start", "end", "strand", "origin_seq"])
    assert filecmp.cmp("tmp", os.path.join("tests", "data", "cds.csv"))
    os.remove('tmp')
    # Export only two columns
    export.export_csv(info, "tmp", ["start", "end"])
    assert filecmp.cmp("tmp", os.path.join("tests", "data", "light_cds.csv"))
    os.remove('tmp')


def test_export_fasta():
    """Test export_fasta function"""
    seq = [
        SeqRecord(
            Seq("MKEIPFDDFINNISKAETPFAVKIDNMLLQGNCLREIKEGASGYSVRMKRHESISW*"),
            id="scaffold_117_243_1",
            description="# origin_seq: scaffold_117 # strand: reverse # start: 228271 # end: 228441"),
        SeqRecord(
            Seq("MKEIPFDDFINNISKAETPFAVKIDNMLLQGNCLREIKEGASGYSVRMKRHESISWOLPDYKSNKGVFKSELISVLIFQV*"),
            id="scaffold_117_243_2",
            description="# origin_seq: scaffold_117 # strand: reverse # start: 228199 # end: 228441"),
        SeqRecord(
            Seq("VANEEKKDFNAMLRKNTDMPKTQIVTDESIIKRYGGERMFFAPPLAYDELMKKVPHGKVVTAEKIREYLAEKNCADFTDPMTAGLFISIAAWASHQREEDITPYWRTLKTDGELNAKYPGGIEAQKKMLEEEGHVIIQKGRKNIRFFVKDYENVLFDLH*"),
            id="scaffold_1220_19_1",
            description="# origin_seq: scaffold_1220 # strand: forward # start: 15451 # end: 15930"),
        SeqRecord(
            Seq("VANEEKKDFNAMLRKNTDMPKTQIVTDESIIKRYGGERMFFAPPLAYDELMKKVPHGKVVTAEKIREYLAEKNCADFTDPMTAGLFISIAAWASHQREEDITPYWRTLKTDGELNAKYPGGIEAQKKMLEEEGHVIIQKGRKNIRFFVKDYENVLFDLHOSDLVKLFL*"),
            id="scaffold_1220_19_2",
            description="# origin_seq: scaffold_1220 # strand: forward # start: 15451 # end: 15957")]
    export.export_fasta(seq, "tmp")
    assert filecmp.cmp("tmp", os.path.join("tests", "data", "pot_pyl_cds.fasta"))
