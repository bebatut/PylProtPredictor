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
    start = int(split_description[1].replace(" ", ""))
    end = int(split_description[2].replace(" ", ""))
    strand = transform_strand(split_description[3].replace(" ", ""))
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
    return (end + 3) < origin_seq_size and (end + 3) < next_cds_limit


def find_stop_codon_pos_in_seq(seq):
    """Find position of STOP codon in a sequence

    :param seq: string sequence of amino acids

    :return: list of position for possible STOP codons in a sequence
    """
    stop_codon_pos = []
    for i in range(len(seq) - 1):
        if seq.startswith('*', i):
            stop_codon_pos.append(i)
    return stop_codon_pos


def translate(seq):
    """Translate a sequence into amino acids while replacing any possible STOP
    codon encoded by TAG by a Pyl amino acid

    :param seq: a Seq object

    :return: string with the corresponding amino acid sequence with the TAG encoded STOP are replaced by Pyl amino acid
    """
    translated_seq = seq.translate()
    str_seq = str(seq)
    n = 3
    codons = [str_seq[i:(i + n)] for i in range(0, len(str_seq), n)]
    for i in find_stop_codon_pos_in_seq(str(translated_seq)):
        if codons[i] != "TAG":
            raise ValueError("Stop codon (different from TAG) found inside the sequence")
        mutable_seq = translated_seq.tomutable()
        mutable_seq[i] = "O"
        translated_seq = mutable_seq.toseq()
    return translated_seq


class CDS:
    'Class to describe a CDS'

    def __init__(
            self, seq_id="", origin_seq=None, origin_seq_id="", start=-1,
            end=-1, strand="forward", seq=None, order=-1, next_cds_limit=-1,
            alternative_ends=[], alternative_cds=[], evalue=10, 
            conserved_cds=None, rejected_cds=[]):
        """Initiate a CDS instance"""
        self.id = seq_id
        self.origin_seq = origin_seq
        self.origin_seq_id = origin_seq_id
        self.start = start
        self.end = end
        self.strand = strand
        self.seq = seq
        self.order = order
        self.next_cds_limit = next_cds_limit
        self.alternative_ends = alternative_ends
        self.alternative_cds = alternative_cds
        self.evalue = evalue
        self.conserved_cds = conserved_cds
        self.rejected_cds = rejected_cds

    def init_from_record(self, record):
        """Initiate a CDS instance with a SeqRecord object

        :param record:
        """
        seq_id, origin_seq_id, start, end, strand = extract_seq_desc(record.description)
        self.set_id(seq_id)
        self.set_origin_seq_id(origin_seq_id)
        self.set_start(start)
        self.set_end(end)
        self.set_strand(strand)
        self.set_seq(record.seq)

    def init_from_dict(self, in_dict):
        """Initiate a CDS instance with a dictionary

        :param in_dict: dictionary with attribute for a CDS object
        """
        self.set_id(in_dict["id"])
        self.set_origin_seq_id(in_dict["origin_seq_id"])
        self.set_start(in_dict["start"])
        self.set_end(in_dict["end"])
        self.set_strand(in_dict["strand"])
        self.set_seq(Seq(in_dict["seq"]))
        self.set_order(in_dict["order"])
        self.set_next_cds_limit(in_dict["next_cds_limit"])
        self.set_alternative_ends(in_dict["alternative_ends"])

        self.reset_alternative_cds()
        for cds_id in in_dict["alternative_cds"]:
            print(cds_id)
            new_cds = CDS()
            new_cds.init_from_dict(in_dict["alternative_cds"][cds_id])
            self.add_alternative_cds(new_cds)

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

    def get_evalue(self):
        """Return the evalue

        :return: evalue of the CDS object
        """
        return self.evalue

    def get_conserved_cds(self):
        """Return the CDS object of the conserved CDS as correct CDS (start, end, sequence)

        :return: CDS object of the conserved CDS
        """
        return self.conserved_cds

    def get_rejected_cds(self):
        """Return a list of the rejected CDS objects as correct CDS (start, end, sequence)

        :return: list of CDS objects
        """
        return self.rejected_cds

    def get_origin_seq_size(self):
        """Return the length of the origin sequence

        :return: int corresponding to the length of the origin sequence
        """
        if not self.has_origin_seq():
            raise ValueError("No origin sequence provided")
        return len(self.get_origin_seq().seq)

    def get_origin_seq_string(self):
        """Return the string of the origin sequence

        :return: string corresponding to the origin sequence
        """
        if not self.has_origin_seq():
            raise ValueError("No origin sequence provided")
        return str(self.get_origin_seq().seq)

    def get_alternative_start(self):
        """Return the list of alternative CDS start

        :return: list of the start of the alternative CDS
        """
        alt_starts = []
        for alt_cds in self.get_alternative_cds():
            alt_starts.append(alt_cds.get_start())
        return alt_starts

    def get_alternative_end(self):
        """Return the list of alternative CDS end

        :return: list of the end of the alternative CDS
        """
        alt_ends = []
        for alt_cds in self.get_alternative_cds():
            alt_ends.append(alt_cds.get_end())
        return alt_ends

    def get_translated_seq(self):
        """Return the translated sequence of the CDS

        :return: SeqRecord object corresponding to the translated sequence
        """
        seq = SeqRecord(
            translate(self.get_seq()),
            id=self.get_id(),
            description=self.export_description())
        return seq

    def get_translated_alternative_seq(self):
        """Return a list of the translated sequences of the alternative sequences

        :return: list of SeqRecord objects
        """
        transl_alt_seq = []
        for alt_cds in self.get_alternative_cds():
            transl_alt_seq.append(alt_cds.get_translated_seq())
        return transl_alt_seq

    def get_seqrecord(self):
        """Return a SeqRecord of the CDS

        :return: SeqRecord
        """
        seq = SeqRecord(
            self.get_seq(),
            id=self.get_id(),
            description=self.export_description())
        return seq 


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
        if strand != "reverse" and strand != "forward":
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

    def set_evalue(self, evalue):
        """Change the evalue

        :param evalue: new evalue
        """
        self.evalue = evalue

    def set_conserved_cds(self, conserved_cds):
        """Change the conserved CDS

        :param conserved_cds: CDS object of the conserved CDS
        """
        self.conserved_cds = conserved_cds

    def reset_alternative_cds(self):
        """Reset the list of alternative cds"""
        self.alternative_cds = []

    def add_alternative_cds(self, alternative_cds):
        """Add an alternative CDS to the list of possible alternative CDS

        :param alternative_cds: a CDS object
        """
        self.alternative_cds.append(alternative_cds)

    def reset_rejected_cds(self):
        """Reset the list of rejected cds"""
        self.rejected_cds = []

    def add_rejected_cds(self, rejected_cds):
        """Add a rejected CDS to the list of rejected CDS

        :param rejected_cds: a CDS object
        """
        self.rejected_cds.append(rejected_cds)

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

    def has_next_cds_limit(self):
        """Test if the limit of next CDS is set

        :return: boolean
        """
        return self.get_next_cds_limit() != -1

    def has_alternative_ends(self):
        """Test if the list of alternative ends is not empty

        :return: boolean
        """
        return len(self.get_alternative_ends()) > 0

    def is_potential_pyl_cds(self):
        """Test if the CDS is a potential PYL CDS: having alternative cds

        :return: boolean
        """
        return len(self.get_alternative_cds()) > 0

    def has_origin_seq(self):
        """Test if the CDS has a origin seq

        :return: boolean
        """
        return self.get_origin_seq() is not None

    def find_next_cds_limit(self, ordered_pred_cds, pred_cds):
        """Determine the end of the next CDS on the strand on the origin sequence

        If the stand is reverse, we need to take the start of the previous CDS

        :param ordered_pred_cds: ordered ids of the CDS on the same strand on the origin sequence
        :param pred_cds: a dictionary with the predicted CDS represented as CDS objects
        """
        if ordered_pred_cds[self.get_order()] != self.get_id():
            raise ValueError("Incorrect order for %s" % (self.get_id()))

        if not self.has_origin_seq():
            raise ValueError("No origin sequence provided")

        origin_seq_size = self.get_origin_seq_size()

        if self.is_reverse_strand():
            next_id = self.get_order() - 1
            next_cds_limit = 0
            if next_id >= 0:
                next_cds_id = ordered_pred_cds[next_id]
                next_cds_limit = pred_cds[next_cds_id].get_start()
            next_cds_limit = (origin_seq_size - next_cds_limit + 1)
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
        if not self.has_origin_seq():
            raise ValueError("No origin sequence provided")

        if not self.has_next_cds_limit():
            raise ValueError("No next CDS limit provided")

        origin_seq = self.get_origin_seq_string()
        origin_seq_size = self.get_origin_seq_size()

        new_end = self.get_end()
        if self.is_reverse_strand():
            new_end = (origin_seq_size - self.get_start() + 1)

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

    def extract_possible_alternative_seq(self):
        """
        Extract the start, end and sequence of different possible sequences for a CDS identified as
        potential PYL CDS
        """
        if not self.has_alternative_ends():
            return

        if not self.has_origin_seq():
            raise ValueError("No origin sequence provided")

        origin_seq = self.get_origin_seq_string()
        origin_seq_size = self.get_origin_seq_size()

        start = self.get_start()
        end = self.get_end()
        seq_id = self.get_id()

        self.reset_alternative_cds()
        count = 1
        for alt_end in self.get_alternative_ends():
            if self.is_reverse_strand():
                new_start = origin_seq_size - alt_end + 1
                new_end = end
                rev_start = origin_seq_size - end + 1
                new_seq = origin_seq[(rev_start - 1):alt_end]
            else:
                new_start = start
                new_end = alt_end
                new_seq = origin_seq[(start - 1):new_end]
            new_cds = CDS(
                seq_id="%s-%s" % (seq_id, count),
                start=new_start,
                end=new_end,
                strand=self.get_strand(),
                origin_seq_id=self.get_origin_seq_id(),
                seq=Seq(new_seq))
            count += 1
            self.add_alternative_cds(new_cds)

    def add_evalue(self, seq_id, evalue):
        """Add evalue to the correct CDS object

        :param seq_id: id of the CDS
        :param evalue: evalue to add
        :param origin_cds_obj: original CDS object
        """
        if seq_id == self.get_id():
            self.set_evalue(evalue)
        else:
            for alt_cds in self.get_alternative_cds():
                if alt_cds.get_id() == seq_id:
                    alt_cds.set_evalue(evalue)

    def identify_cons_rej_cds(self):
        """
        Identify which alternative CDS has conversed or rejected based on the 
        evalue
        """
        ref_evalue = self.get_evalue()
        self.reset_rejected_cds()
        self.set_conserved_cds(self)
        for alt_cds in self.get_alternative_cds():
            evalue = alt_cds.get_evalue()
            if evalue < ref_evalue:
                ref_evalue = evalue
                self.add_rejected_cds(self.get_conserved_cds())
                self.set_conserved_cds(alt_cds)
            else:
                self.add_rejected_cds(alt_cds)

    def export_description(self):
        """Export the description of the CDS

        :return: string with the description
        """
        desc = "# origin_seq: %s # strand: %s # start: %s # end: %s" % (
            self.get_origin_seq_id(),
            self.get_strand(),
            self.get_start(),
            self.get_end())
        return desc

    def export_to_dict(self):
        """Export the object to CDS

        :return: dict corresponding to CDS object
        """
        cds_id = self.get_id()
        d = {cds_id : {
            'id': self.get_id(),
            'origin_seq_id': self.get_origin_seq_id(),
            'start': self.get_start(),
            'end': self.get_end(),
            'strand': self.get_strand(),
            'seq': str(self.get_seq()),
            'order': self.get_order(),
            'next_cds_limit': self.get_next_cds_limit(),
            'alternative_ends': self.get_alternative_ends(),
            'alternative_cds': {}}
        }    
        
        if self.has_alternative_ends():
            for alt_cds in self.get_alternative_cds():
                d[cds_id]['alternative_cds'].update(alt_cds.export_to_dict())

        return d

