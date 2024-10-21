from io import StringIO
from pandas import DataFrame
import pytest
import textwrap
from commec.tools.blast_tools import _split_by_tax_id, readblast

@pytest.fixture
def blast_df():
    """
    Return a dataframe containing 3 BLAST hits, 2 with multiple taxids, 1 of which is invalid and 1
    of which is a synthetic taxid
    """
    blast_to_parse = textwrap.dedent(
        """\
        # BLASTX 2.15.0+
        # Query: NC_TEST
        # Database: /root/commec-dbs/mock
        #query acc.	subject title	subject acc.	subject tax ids	evalue	bit score	% identity	query length	q. start	q. end	subject length	s. start	s. end
        # 3 hits found
        BT_01	SUBJECT	SUBJECT_ACC	2371;644357	0.0	BITSCORE	99.999	300	101	200	500	1	100
        BT_01	SUBJECT	SUBJECT_ACC	10760;110011001100	0.0	BITSCORE	99.999	300	25	80	500	1	100
        BT_01	SUBJECT	SUBJECT_ACC	32630	0.0	BITSCORE	99.999	300	275	300	500	1	100
        """
    )
    return readblast(StringIO(blast_to_parse))

def test_split_by_tax_id(blast_df: DataFrame):
    assert len(blast_df) == 3
    split_blast = _split_by_tax_id(blast_df)
    assert len(split_blast) == 5
    expected_tax_ids = {2371, 644357, 10760, 110011001100, 32630}
    assert set(split_blast["subject tax ids"]) == expected_tax_ids
