#!/usr/bin/env python

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Data import CodonTable


def extract_seq_desc(desc):
    """Extract from description the seq id, the origin sequence, start, end and
    strand from a predicted CDS

    :param desc: description of a prediced CDS with Prodigal

    :return: id of predicted CDS
    :return: id of the origin sequence
    :return: start position of the predicted CDS
    :return: end position of the predicted CDS
    :return: strand of the predicted CDS
    """
    split_description = desc.split("#")
    seq_id = split_description[0][:-1]
    origin_seq = "_".join(seq_id.split("_")[:-1])
    if seq_id.find("|") != -1:
        seq_id = "cds_%s" % (seq_id.split("_")[-1])
    start = int(split_description[1].replace(" ",""))
    end = int(split_description[2].replace(" ",""))
    strand = transform_strand(split_description[3].replace(" ",""))
    return seq_id, origin_seq, start, end, strand


def transform_strand(strand_id):
    """Transform strand from numerical value to string value

    :param strand_id: numerical value to represent a strand (1 or -1)

    :return: string value (forward or reverse) for the strand
    """
    if strand_id == "-1":
        return 'reverse'
    elif strand_id == "1":
        return 'forward'
    else:
        raise ValueError("Wrong strand_id: {}".format(strand_id))


def test_to_continue(end, origin_seq_size, next_cds_limit):
    """Test if possible to extract next codon: position still in the genome and
    lower than the end of the next CDS

    :param end: int corresponding to the current end
    :param origin_seq_size: size of the origin sequence
    :param next_cds_limit: int corresponding to the end of the next CDS on the same strand

    :return: boolean
    """
    return (end+3) < origin_seq_size and (end+3) < next_cds_limit


class CDS:
    'Class to describe a CDS'

    def __init__(self):
        """Initiate a CDS instance"""
        self.id = ""
        self.origin_seq = None
        self.origin_seq_id = ""
        self.start = -1
        self.end = -1
        self.strand = "forward"
        self.seq = None
        self.order = -1
        self.next_cds_limit = -1
        self.alternative_ends = []
        self.alternative_cds = []


    def init_from_record(record):
        """Initiate a CDS instance with a SeqRecord object"""
        seq_id, origin_seq_id, start, end, strand = extract_seq_desc(record.description)
        self.set_id(seq_id)
        self.set_origin_seq_id(origin_seq_id)
        self.set_start(start)
        self.set_end(end)
        self.set_strand(strand)
        self.set_seq(record.seq)


    def get_id(self):
        """Return the id of the CDS

        :return: string corresponding to the id
        """
        return self.id


    def get_origin_seq_id(self):
        """Return the id of origin seq of the CDS

        :return: string corresponding to the origin seq
        """
        return self.origin_seq_id


    def get_origin_seq(self):
        """Return the SeqRecord object corresponding to the origin seq of the CDS

        :return: SeqRecord object
        """
        return self.origin_seq


    def get_start(self):
        """Return the start position of the CDS on the origin sequence

        :return: int corresponding to the start position
        """
        return self.start


    def get_end(self):
        """Return the end position of the CDS on the origin sequence

        :return: int corresponding to the end position
        """
        return self.end


    def get_strand(self):
        """Return the strand of the CDS on the origin sequence

        :return: string corresponding to the strand (forward or reverse)
        """
        return self.strand


    def get_seq(self):
        """Return the sequence of the CDS

        :return: string with the sequence
        """
        return self.seq


    def get_order(self):
        """Return the order of the CDS on the strand on the origin sequence

        :return: int corresponding to the order value
        """
        return self.order


    def get_next_cds_limit(self):
        """Return the end or start (if reverse strand) of next CDS on the strand on the origin sequence

        :return: int corresponding to the end or start of the next CDS
        """
        return self.next_cds_limit


    def get_alternative_ends(self):
        """Return the list of possible alternative ends if the CDS is ending with TAG STOP codon

        :return: list of int corresponding to the alternative ends
        """
        return self.alternative_ends


    def get_alternative_cds(self):
        """Return the list of possible alternative CDS if the CDS is ending with TAG STOP codon

        :return: list of CDS object
        """
        return self.alternative_cds


    def set_id(self, seq_id):
        """Change the id of the CDS

        :param order: new id value
        """
        self.id = seq_id


    def set_origin_seq_id(self, origin_seq_id):
        """Change the id of the origin sequence of the CDS

        :param order: new origin seq id value
        """
        self.origin_seq_id = origin_seq_id


    def set_start(self, start):
        """Change the start position of the CDS

        :param start: new start value (int)
        """
        self.start = start


    def set_end(self, end):
        """Change the end position of the CDS

        :param end: new end value (int)
        """
        self.end = end


    def set_strand(self, strand):
        """Change the strand value of the CDS

        :param end: new strand (forward or reverse)
        """
        if strand != "reverse" or strand != "forward":
            raise ValueError("Incorrect strand value: %s" % (strand))
        self.strand = strand


    def set_seq(self, seq):
        """Change the sequence object of the CDS

        :param seq: new Seq object with the sequence of the CDS
        """
        self.seq = seq


    def set_order(self, order):
        """Change the order of the CDS on the strand on the origin sequence

        :param order: new order value
        """
        self.order = order


    def set_next_cds_limit(self, next_cds_limit):
        """Change the end or start (if reverse strand) of next CDS on the strand on the origin sequence

        :param next_cds_limit: int corresponding to the end or start of next CDS
        """
        self.next_cds_limit = next_cds_limit


    def set_origin_seq(self, origin_seq):
        """Change the SeqRecord object corresponding to the origin seq of the CDS

        :param origin_seq: SeqRecord object
        """
        if self.get_strand() == "reverse":
            self.origin_seq = origin_seq.reverse_complement()
            self.origin_seq.description = origin_seq.description
            self.origin_seq.id = origin_seq.id
            self.origin_seq.name = origin_seq.name
        else:
            self.origin_seq = origin_seq


    def set_alternative_ends(self, alternative_ends):
        """Change the list of alternative ends

        :param alternative_ends: list of int corresponding to the new alternative ends
        """
        self.alternative_ends = alternative_ends


    def add_alternative_cds(self, alternative_cds):
        """Add an alternative CDS to the list of possible alternative CDS

        :param alternative_cds: a CDS object
        """
        self.alternative_cds.apend(alternative_cds)


    def is_reverse_strand(self):
        """Test if the strand is reverse

        :return: boolean
        """
        return self.get_strand() == "reverse"


    def is_tag_ending_seq(self):
        """Test if the sequence is ending with TAG STOP codon 

        :return: boolean
        """
        return self.get_seq().endswith("TAG")


    def has_alternative_ends(self):
        """Test if the list of alternative ends is not empty

        :return: boolean
        """
        return len(self.get_alternative_ends()) > 0


    def find_next_cds_limit(self, ordered_pred_cds, pred_cds):
        """Determine the end of the next CDS on the strand on the origin sequence

        If the stand is reverse, we need to take the start of the previous CDS
        
        :param ordered_pred_cds: ordered ids of the CDS on the same strand on the origin sequence
        :param pred_cds: a dictionary with the predicted CDS represented as CDS objects 
        """
        if ordered_pred_cds[self.get_order()] != self.get_id():
            raise ValueError("Incorrect order for %s" % (self.get_id()))

        if self.get_origin_seq() is None:
            raise ValueError("No origin sequence provided")

        origin_seq_size = len(self.get_origin_seq().seq)
        
        if self.is_reverse_strand():
            next_id = self.get_order() - 1
            next_cds_limit = 0
            if next_id >= 0:
                next_cds_id = ordered_pred_cds[next_id]
                next_cds_limit = pred_cds[next_cds_id].get_start()
            next_cds_limit = (origin_seq_size-next_cds_limit+1)
        else:
            next_id = self.get_order() + 1
            next_cds_limit = origin_seq_size
            if next_id < len(ordered_pred_cds):
                next_cds_id = ordered_pred_cds[next_id]
                next_cds_limit = pred_cds[next_cds_id].get_end()

        self.set_next_cds_limit(next_cds_limit)


    def find_alternative_ends(self):
        """
        Find alternative ends (on the same ORF) for a CDS until the next found STOP 
        codon on the genome (or its complement if the CDS is on the reverse strand)
        """
        if self.get_origin_seq() is None:
            raise ValueError("No origin sequence provided")

        origin_seq = str(self.get_origin_seq().seq)
        origin_seq_size = len(origin_seq.seq)

        if self.get_next_cds_limit == -1:
            raise ValueError("No next CDS limit provided")

        new_end = self.get_end()
        if self.is_reverse_strand():
            new_end = (origin_seq_size-self.get_start()+1)

        stop_codons = CodonTable.unambiguous_dna_by_id[1].stop_codons
        to_continue = test_to_continue(new_end, origin_seq_size, self.get_next_cds_limit())
        new_ends = []
        while to_continue:
            codon = origin_seq[new_end:(new_end + 3)]
            new_end += 3
            if codon not in stop_codons:
                to_continue = test_to_continue(new_end, origin_seq_size, self.get_next_cds_limit())
            else:
                new_ends.append(new_end)
                if codon != 'TAG':
                    to_continue = False
                else:
                    to_continue = test_to_continue(new_end, origin_seq_size, self.get_next_cds_limit())
        self.set_alternative_ends(new_ends)


    def extract_possible_seq(self):
        """
        Extract the start, end and sequence of different possible sequences for a CDS identified as
        potential PYL CDS
        """
        if not self.has_alternative_ends():
            return
        
        if self.get_origin_seq() is None:
            raise ValueError("No origin sequence provided")

        origin_seq = str(self.get_origin_seq().seq)
        origin_seq_size = len(origin_seq.seq)

        for new_end in self.get_alternative_ends():
            new_start = self.get_start()
            new_seq = origin_seq[(start-1):new_end]
            if strand == "reverse":
                new_start = origin_seq_size - new_end + 1
                new_end = origin_seq_size - start + 1
            new_cds = CDS()
            new_cds.set_start(new_start)
            new_cds.set_end(new_end)
            new_cds.set_strand(self.get_strand())
            new_cds.set_origin_seq_id(self.get_origin_seq_id())
            new_cds.set_seq(new_seq)

