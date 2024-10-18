from io import StringIO
import os
import pandas as pd
import pytest
from Bio import SeqIO
from commec.screeners.fetch_nc_bits import (
    get_ranges_with_no_hits,
    write_nc_sequences,
    fetch_noncoding_regions,
)


@pytest.mark.parametrize(
    "hits, nc_ranges",
    [
        # Two protein hits, no noncoding regions > 50bp
        ([(1, 50), (100, 150), (175, 299)], []),
        # One protein hit, < 50bp nocoding regions on the ends
        ([(50, 251)], []),
        # One protein hit, > 50bp nocoding regions on the ends
        ([(51, 250)], [[1, 50], [251, 300]]),
        # Three protein hits, one noncoding region >50bp
        (
            [(1, 40), (140, 265), (300, 349)],
            [[41, 139]],
        ),
    ],
)
def test_get_ranges_with_no_hits(hits, nc_ranges):
    """
    Test the BLAST hits are successfully converted into noncoding ranges.
    """
    def _create_mock_blast_df_from(hits):
        data = {
            "q. start": [hit[0] for hit in hits],
            "q. end": [hit[1] for hit in hits],
            "query length": [300] * len(hits),
        }
        df = pd.DataFrame(data)
        return df.reset_index(drop=True)  # This adds a numeric index

    blast_df = _create_mock_blast_df_from(hits)
    assert get_ranges_with_no_hits(blast_df) == nc_ranges


def test_write_nc_sequences(tmp_path):
    """
    Test that sequences are successfully written to an output file based on noncoding regions.
    """
    desc = "TEST_FETCH_NC"
    seq = "".join([str(i % 10) for i in range(300)])
    record = list(SeqIO.parse(StringIO(f">{desc}\n{seq}\n"), "fasta"))[0]
    nc_ranges = [[5, 59], [91, 170]]
    outfile = tmp_path / "output.fasta"

    write_nc_sequences(nc_ranges, record, outfile)

    content = outfile.read_text()
    assert f">{desc} 5-59\n{seq[4:59]}\n" in content
    assert f">{desc} 91-170\n{seq[90:170]}\n" in content


INPUT_FASTA = os.path.join(os.path.dirname(__file__), "test_data/fetch_nc_input.fasta")
INPUT_BLAST = os.path.join(os.path.dirname(__file__), "test_data/fetch_nc_input.blastx")
OUTPUT_FILE = os.path.join(
    os.path.dirname(__file__), "test_data/fetch_nc_input.blastx.noncoding.fasta"
)
EXPECTED_FILE = os.path.join(
    os.path.dirname(__file__), "test_data/fetch_nc_expected.fasta"
)
def test_main():
    """
    Functional run of entire fetch_nc_bits, on a mock blast, and mock input fasta.
    """
    fetch_noncoding_regions(INPUT_BLAST, INPUT_FASTA)

    with open(OUTPUT_FILE, "r", encoding="utf-8") as file1, open(
        EXPECTED_FILE, "r", encoding="utf-8"
    ) as file2:
        for line1, line2 in zip(file1, file2):
            # Strip the lines to ignore trailing spaces or newlines
            assert line1.strip() == line2.strip()
