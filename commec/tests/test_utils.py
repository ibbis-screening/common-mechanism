from commec.utils import is_protein_specific, is_mostly_standard_bases

def test_is_protein_specific():
    # has IUPAC characters unique to a protein
    sequence = 'ADEYFMMGDSIVLKIVLHGGDPYFEPRMHNITS'
    assert is_protein_specific(sequence) is True

    # nucleotide record (could also be an unusual protein)
    sequence = 'gtagacaacaaattcaacaaagaacaacaacgt'
    assert is_protein_specific(sequence) is False

def test_is_mostly_standard_bases():
    # has IUPAC characters unique to a protein
    sequence = 'ADEYFMMGDSIVLKIVLHGGDPYFEPRMHNITS'
    assert is_mostly_standard_bases(sequence) is False

    # nucleotide record (could also be an unusual protein)
    sequence = 'gtagacaacaaattcaacaaagaacaacaacgt'
    assert is_mostly_standard_bases(sequence) is True