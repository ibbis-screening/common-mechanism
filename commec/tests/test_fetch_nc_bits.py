from unittest.mock import call, patch, mock_open
import pytest
import pandas as pd
from io import StringIO
from Bio import SeqIO
from commec.fetch_nc_bits import get_ranges_with_no_hits, write_nc_sequences


# Mock sequence using digits 0-9
MOCK_SEQ_LENGTH = 300
MOCK_SEQ = "".join([str(i % 10) for i in range(MOCK_SEQ_LENGTH)])
MOCK_FASTA = f">test_sequence\n{MOCK_SEQ}\n"


@pytest.fixture
def mock_fasta_file():
    return StringIO(MOCK_FASTA)


def create_mock_blast_df(hits):
    data = {
        "q. start": [hit[0] for hit in hits],
        "q. end": [hit[1] for hit in hits],
        "query length": [300] * len(hits),
    }
    df = pd.DataFrame(data)
    return df.reset_index(drop=True)  # This adds a numeric index


@pytest.mark.parametrize(
    "hits, nc_ranges",
    [
        ([(1, 50), (100, 150), (175, 299)], []),  # Two protein hits, no gaps > 50bp
        ([(50, 251)], []),  # One protein hit, < 50bp nc on the ends
        ([(51, 250)], [[1, 50], [251, 300]]),  # One protein hit, > 50bp nc on the ends
        (
            [(1, 40), (140, 265), (300, 349)],
            [[41, 139]],
        ),  # Three protein hits, one gap >50bp
    ],
)
def test_get_ranges_with_no_hits(hits, nc_ranges):
    blast_df = create_mock_blast_df(hits)
    assert get_ranges_with_no_hits(blast_df) == nc_ranges


def test_write_nc_sequences_partial_sequence(mock_fasta_file):
    nc_ranges = [[5, 59], [91, 170]]
    record = list(SeqIO.parse(mock_fasta_file, "fasta"))[0]

    with patch("commec.fetch_nc_bits.open", new_callable=mock_open) as mock_file:
        # Set up the mock to return our StringIO object when opened for reading
        mock_file.return_value.__enter__.return_value.read.return_value = MOCK_FASTA

        write_nc_sequences(nc_ranges, record, "output.fasta")

        mock_file.assert_any_call("output.fasta", "w", encoding="utf-8")
        write_calls = mock_file.return_value.__enter__.return_value.writelines.call_args_list
        assert len(write_calls) == 1, f"Expected 1 writelines call, but got {len(write_calls)}"

        expected_call = call(
                f">test_sequence test_sequence 5-59\n{''.join([str(i % 10) for i in range(4, 58)])}\n" +
                f">test_sequence test_sequence 91-170\n{''.join([str(i % 10) for i in range(90, 169)])}\n"
        )

        assert expected_call == write_calls
