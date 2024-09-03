""" 
Unit test for ensuring that the databases are being called without errors.
"""
import os
import pytest
from commec.databases.blastdmnd_db import DiamondDataBase
from commec.databases.blastn_db import BlastNDataBase
from commec.databases.blastx_db import BlastXDataBase
from commec.databases.hmm_db import HMMDataBase
from commec.databases.cmscan_db import CmscanDataBase

INPUT_QUERY = os.path.join(os.path.dirname(__file__),"test_data/single_record.fasta")
DATABASE_DIRECTORY = os.path.join(os.path.dirname(__file__),"test_dbs")

databases_to_implement = [
    [DiamondDataBase,   "nr_dmnd",     "nr.dmnd"],
    [BlastNDataBase,    "nt_blast",    "nt"],
    [BlastXDataBase,    "nr_blast",    "nr"],
    [HMMDataBase,       "benign_db",   "benign.hmm"],
    #[CmscanDataBase,    "biorisk_db",  "biorisk.hmm"],
]

@pytest.mark.parametrize("input_db", databases_to_implement)
def test_database_can_run(input_db):
    """
    Opens a database object on a test database, and runs the test query on it.
    Fails if commec environment is not setup correctly, or if the database object
    defaults are invalid etc.

    Something similar to this would be useful to be run 
    instead of --help during the conda recipe checks.
    """

    db_dir = os.path.join(DATABASE_DIRECTORY, input_db[1])
    db_file = os.path.join(db_dir, input_db[2])

    output_file = "db.out"

    new_db = input_db[0](db_dir, db_file, INPUT_QUERY, output_file)
    new_db.screen()
    assert new_db.check_output()

    if os.path.isfile(output_file):
        os.remove(output_file)
 