import os
from unittest.mock import mock_open, patch
import pytest
import pandas as pd
from Bio import SeqIO
from commec.fetch_nc_bits import get_ranges_with_no_hits, write_nc_sequences


MOCK_SEQ_LENGTH = 300
MOCK_SEQ = "ATGC" * int(MOCK_SEQ_LENGTH / 4)

def create_mock_blast_df(hits):
    data = {
        "q. start": [hit[0] for hit in hits],
        "q. end": [hit[1] for hit in hits],
        "query length": [MOCK_SEQ_LENGTH] * len(hits)
    }
    df = pd.DataFrame(data)
    return df.reset_index(drop=True)  # This adds a numeric index

@pytest.mark.parametrize("hits, nc_ranges", [
    ([(1, 50), (100, 150), (175, 299)], []),  # Two protein hits, no gaps > 50bp
    ([(50, 251)], []),                        # One protein hit, < 50bp nc on the ends
    ([(51, 250)], [[1, 50], [251, 300]]),     # One protein hit, > 50bp nc on the ends
])
def test_get_ranges_with_no_hits(hits, nc_ranges):
    blast_df = create_mock_blast_df(hits)
    assert get_ranges_with_no_hits(blast_df) == nc_ranges
