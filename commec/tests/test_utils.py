from commec.utils import is_protein_specific

def test_is_protein_specific():
    sequence = ''
    assert is_protein_specific(sequence) is False

    # has IUPAC characters unique to a protein
    sequence = 'ADEYFMMGDSIVLKIVLHGGDPYFEPRMHNITS'
    assert is_protein_specific(sequence) is True

    # nucleotide record (could also be an unusual protein)
    sequence = 'ADEYFMMGDSIVLKIVLHGGDPYFEPRMHNITS'
    assert is_protein_specific(sequence) is True
