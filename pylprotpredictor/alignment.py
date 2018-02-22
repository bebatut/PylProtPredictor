class Alignment:
    'Class to describe a DIAMOND alignment'

    def __init__(
            self, sseqid="", pident=0, length=0, mismatch=0, gapopen=0, qstart=0,
            qend=0, sstart=0, send=0, evalue=10, bitscore=0):
        """Initiate a CDS instance"""
        self.sseqid = sseqid
        self.pident = pident
        self.length = length
        self.mismatch = mismatch
        self.gapopen = gapopen
        self.qstart = qstart
        self.qend = qend
        self.sstart = sstart
        self.send = send
        self.evalue = evalue
        self.bitscore = bitscore

    def init_from_search_report_row(self, row):
        """Initiate an Alignment instance with a row extracted from a
        BLAST/DIAMOND table

        :param row: a pandas row
        """
        self.set_sseqid(row[2])
        self.set_pident(float(row[3]))
        self.set_length(int(row[4]))
        self.set_mismatch(int(row[5]))
        self.set_gapopen(int(row[6]))
        self.set_qstart(int(row[7]))
        self.set_qend(int(row[8]))
        self.set_sstart(int(row[9]))
        self.set_send(int(row[10]))
        self.set_evalue(float(row[11]))
        self.set_bitscore(float(row[12]))

    def get_sseqid(self):
        """Return query sequence id

        :return: string
        """
        return self.sseqid

    def get_pident(self):
        """Return percentage of identical matches

        :return: float
        """
        return self.pident

    def get_length(self):
        """Return alignment length

        :return: int
        """
        return self.length

    def get_mismatch(self):
        """Return number of mismatches

        :return: int
        """
        return self.mismatch

    def get_gapopen(self):
        """Return number of gap openings

        :return: int
        """
        return self.gapopen

    def get_qstart(self):
        """Return query sequence id

        :return: query sequence id
        """
        return self.qstart

    def get_qend(self):
        """Return end of alignment in query

        :return: int
        """
        return self.qend

    def get_sstart(self):
        """Return start of alignment in subject

        :return: int
        """
        return self.sstart

    def get_send(self):
        """Return end of alignment in subject

        :return: int
        """
        return self.send

    def get_evalue(self):
        """Return expect value

        :return: float
        """
        return self.evalue

    def get_bitscore(self):
        """Return bit score

        :return: query sequence id
        """
        return self.bitscore

    def set_sseqid(self, sseqid):
        """Modify query Seq-id

        :param sseqid: string
        """
        self.sseqid = sseqid

    def set_pident(self, pident):
        """Modify percentage of identical matches

        :param pident: float
        """
        self.pident = pident

    def set_length(self, length):
        """Modify alignment length

        :param length: int
        """
        self.length = length

    def set_mismatch(self, mismatch):
        """Modify number of mismatches

        :param mismatch: int
        """
        self.mismatch = mismatch

    def set_gapopen(self, gapopen):
        """Modify number of gap openings

        :param gapopen: string
        """
        self.gapopen = gapopen

    def set_qstart(self, qstart):
        """Modify start of alignment in query

        :param qstart: string
        """
        self.qstart = qstart

    def set_qend(self, qend):
        """Modify end of alignment in query

        :param qend: string
        """
        self.qend = qend

    def set_sstart(self, sstart):
        """Modify start of alignment in subject

        :param sstart: int
        """
        self.sstart = sstart

    def set_send(self, send):
        """Modify end of alignment in subject

        :param send: int
        """
        self.send = send

    def set_evalue(self, evalue):
        """Modify evalue

        :param evalue: float
        """
        self.evalue = evalue

    def set_bitscore(self, bitscore):
        """Modify bit score

        :param bitscore: int
        """
        self.bitscore = bitscore
