""" 
Unit test for ensuring that the databases are being called without errors.
Will fail if databases have not been installed as expected, with correct versions.
"""
import os
import pytest
from commec.tools.blastdmnd_handler import DiamondHandler
from commec.tools.blastn_handler import BlastNHandler
from commec.tools.blastx_handler import BlastXHandler
from commec.tools.hmm_handler import HMMHandler
from commec.tools.cmscan_handler import CmscanHandler

INPUT_QUERY = os.path.join(os.path.dirname(__file__),"test_data/single_record.fasta")
DATABASE_DIRECTORY = os.path.join(os.path.dirname(__file__),"test_dbs")

# test each database type.
databases_to_implement = [
    [DiamondHandler,   "nr_dmnd",     "nr.dmnd"],
    [BlastNHandler,    "nt_blast",    "nt"],
    [BlastXHandler,    "nr_blast",    "nr"],
    [HMMHandler,       "benign_db",   "benign.hmm"],
    [CmscanHandler,    "benign_db",   "benign.cmscan"],
]

@pytest.mark.parametrize("input_db", databases_to_implement)
def test_database_can_run(input_db, request):
    """
    Opens a database object on a test database, and runs the test query on it.
    Fails if commec environment is not setup correctly, or if the database object
    defaults are invalid etc.

    Something similar to this would be useful to be run 
    instead of --help during the conda recipe checks.
    """

    create_exemplare : bool = request.config.getoption("--gen-examples")
    db_dir = os.path.join(DATABASE_DIRECTORY, input_db[1])
    db_file = os.path.join(db_dir, input_db[2])
    exemplar_file = db_file+".exemplar.out"
    output_file = "db.out"

    if create_exemplare:
        db = input_db[0](db_file, INPUT_QUERY, exemplar_file)
        db.screen()
        assert db.check_output()

    new_db = input_db[0](db_file, INPUT_QUERY, output_file)
    new_db.screen()
    assert new_db.check_output()

    version : str = new_db.get_version_information()
    assert not version == "Version information retrieval error."

    # Read output, and expected output, and compare line by line
    with (open(output_file, 'r', encoding = "utf-8") as f1, 
          open(exemplar_file, 'r', encoding = "utf-8") as f2):
        for line1, line2 in zip(f1, f2):
            if line1.startswith("# Option settings:"):
                continue
            if line1.startswith("# Query file:"):
                continue
            if line1.startswith("# Target file:"):
                continue
            if line1.startswith("# Current dir:"):
                continue
            if line1.startswith("# Database:"):
                continue
            if line1.startswith("# Date:"):
                continue
            assert line1 == line2
    
    if os.path.isfile(output_file):
        os.remove(output_file)
