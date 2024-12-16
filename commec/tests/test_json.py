import pytest
from commec.config.json_io import *
from commec.tools.search_handler import SearchToolVersion
from dataclasses import asdict

@pytest.fixture
def test_screendata():
    '''Fixture to provide the ScreenData for testing.'''
    return ScreenData(
        #recommendation="PASS",
        commec_info = CommecRunInformation(
            commec_version="0.1.2",
            json_output_version=JSON_COMMEC_FORMAT_VERSION,
            biorisk_database_info=SearchToolVersion("HMM 0.0.0","DB 0.0.0"),
            protein_database_info=SearchToolVersion("Blast 0.0.0","DB 0.0.0"),
            nucleotide_database_info=SearchToolVersion("Blast 0.0.0","DB 0.0.0"),
            benign_protein_database_info=SearchToolVersion("Blast 0.0.0","DB 0.0.0"),
            benign_rna_database_info=SearchToolVersion("Blast 0.0.0","DB 0.0.0"),
            benign_synbio_database_info=SearchToolVersion("Blast 0.0.0","DB 0.0.0"),
            time_taken="00:00:00:00",
            date_run="1.1.2024",
        ),
        queries= [
            QueryData(
                query="Query1",
                length=10,
                sequence="ABCDEFGHIJ",
                recommendation = CommecRecommendationContainer(),
                summary_info = CommecSummaryStatistics(),
                hits = [
                    HitDescription(
                        recommendation=CommecScreenStepRecommendation(CommecRecommendation.WARN, CommecScreenStep.BIORISK),
                        name="ImportantProtein1",
                        annotations = {"domain" : ["Bacteria"]},
                        ranges = [
                            MatchRange(
                                e_value = 0.0,
                                match_start = 0,
                                match_end = 10,
                                query_start = 0,
                                query_end = 10
                            )
                        ]
                    )
                ]
            )
        ],
    )

@pytest.fixture
def empty_screendata():
    '''Fixture to provide the ScreenData for testing.'''
    return ScreenData()

@pytest.mark.parametrize("test_data_fixture",["test_screendata", "empty_screendata"])
def test_json_io(tmp_path, request, test_data_fixture):
    ''' Test to ensure that read/write for JSON ScreenData I/O is working correctly.'''
    test_data = request.getfixturevalue(test_data_fixture)
    json_filename1 = tmp_path / "testread1.json"
    json_filename2 = tmp_path / "testread2.json"
    encode_screen_data_to_json(test_data, json_filename1)
    test_data_retrieved = get_screen_data_from_json(json_filename1)
    encode_screen_data_to_json(test_data_retrieved, json_filename2)
    test_data_retrieved_twice = get_screen_data_from_json(json_filename2)

    # Convert both original and retrieved data to dictionaries and compare
    assert asdict(test_data) == asdict(test_data_retrieved), (
        f"JSON Write/Read interpreter failed.\n"
        f"Test JSON Reference data: \n{asdict(test_data)}\n"
        f"Test JSON output data: \n{asdict(test_data_retrieved)}"
    )

    # Convert both original and retrieved data to dictionaries and compare
    assert asdict(test_data) == asdict(test_data_retrieved_twice), (
        f"JSON Write/Read/Write/Read interpreter failed.\n"
        f"Test JSON Reference data: \n{asdict(test_data)}\n"
        f"Test JSON output data: \n{asdict(test_data_retrieved)}"
    )

def test_erroneous_info(tmp_path, test_screendata):
    ''' Test to ensure that read/write for JSON ScreenData I/O is working correctly.'''
    test_data = test_screendata
    json_filename3 = tmp_path / "testread3.json"
    json_filename4 = tmp_path / "testread4.json"

    encode_screen_data_to_json(test_data, json_filename3)
    test_data_retrieved = get_screen_data_from_json(json_filename3)

    # Add erroneous information
    test_data_dict = asdict(test_data_retrieved)
    test_data_dict["ExtraStuff1"] = "ExtraBitStuff1"
    test_data_dict["queries"][0]["ExtraStuff2"] = "ExtraBitStuff2"
    test_data_dict["queries"][0]["hits"][0]["ranges"].append("ExtraStuff3")
    test_data_dict["queries"][0]["hits"][0]["ranges"].append({"ExtraDictStuff4" : 9999})
    test_data_dict2 = encode_dict_to_screen_data(test_data_dict)
    encode_screen_data_to_json(test_data_dict2, json_filename4)
    test_data_retrieved = get_screen_data_from_json(json_filename4)

    # Convert both original and retrieved data to dictionaries and compare
    assert asdict(test_data) == asdict(test_data_retrieved), (
        f"JSON Write/Read interpreter failed.\n"
        f"Test JSON Reference data: \n{asdict(test_data)}\n\n\n\n"
        f"Test JSON output data: \n{asdict(test_data_retrieved)}\n\n\n\n"
    )

def test_recommendation_ordering():
    assert CommecRecommendation.PASS.importance < CommecRecommendation.FLAG.importance
    assert compare(CommecRecommendation.PASS, CommecRecommendation.FLAG) == CommecRecommendation.FLAG

def test_adding_data_to_existing():
    """
    Tests to ensure the mutability of writing to queries is working as expected.
    """
    def write_info(input_query : QueryData):
        input_query.recommendation.biorisk_screen = CommecRecommendation.PASS
    
    new_screen_data = ScreenData()
    new_screen_data.queries.append(QueryData("test01", 10, "ATGCATGCAT", CommecRecommendation.FLAG))
    write_query = new_screen_data.get_query("test01")
    write_info(write_query)
    assert new_screen_data.queries[0].recommendation.biorisk_screen == CommecRecommendation.PASS
