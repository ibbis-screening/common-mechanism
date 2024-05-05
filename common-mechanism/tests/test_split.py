import pytest
from split import clean_description

@pytest.mark.parametrize("description, expected", [
    ("BBa_K620001_P_22737_Coding_\"WT-F87A_(p450)\"", "BBa_K620001_P_22737_Coding_WT-F87A_p450"),
    ("long description" * 20, "longdescription" * 10),
])
def test_clean_description(description, expected):
    assert clean_description(description) == expected
