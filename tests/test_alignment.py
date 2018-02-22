from pylprotpredictor import alignment


al = alignment.Alignment("sseqid", 99, 100, 1, 1, 10, 200, 30, 220, 0.0005, 333.1)


def test_get_sseqid():
    assert al.get_sseqid() == "sseqid"


def test_get_pident():
    assert al.get_pident() == 99


def test_get_length():
    assert al.get_length() == 100


def test_get_mismatch():
    assert al.get_mismatch() == 1


def test_get_gapopen():
    assert al.get_gapopen() == 1


def test_get_qstart():
    assert al.get_qstart() == 10


def test_get_qend():
    assert al.get_qend() == 200


def test_get_sstart():
    assert al.get_sstart() == 30


def test_get_send():
    assert al.get_send() == 220


def test_get_evalue():
    assert al.get_evalue() == 0.0005


def test_get_bitscore():
    assert al.get_bitscore() == 333.1


def test_set_sseqid():
    val = "tt"
    al.set_sseqid(val)
    assert al.get_sseqid() == val


def test_set_pident():
    val = 99.2
    al.set_pident(val)
    assert al.get_pident() == val


def test_set_length():
    val = 2000
    al.set_length(val)
    assert al.get_length() == val


def test_set_mismatch():
    val = 3
    al.set_mismatch(val)
    assert al.get_mismatch() == val


def test_set_gapopen():
    val = 2
    al.set_gapopen(val)
    assert al.get_gapopen() == val


def test_set_qstart():
    val = 30
    al.set_qstart(val)
    assert al.get_qstart() == val


def test_set_qend():
    val = 333
    al.set_qend(val)
    assert al.get_qend() == val


def test_set_sstart():
    val = 5500
    al.set_sstart(val)
    assert al.get_sstart() == val


def test_set_send():
    val = 5550
    al.set_send(val)
    assert al.get_send() == val


def test_set_evalue():
    val = 0.01
    al.set_evalue(val)
    assert al.get_evalue() == val


def test_set_bitscore():
    val = 3333.2
    al.set_bitscore(val)
    assert al.get_bitscore() == val
