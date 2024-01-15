import os
from setuptools import setup
import setuptools

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name = "CommonMechanism",
    version = "0.0.1",
    author = "Nicole Wheeler",
    author_email = "N.Wheeler@bham.ac.uk",
    description = ("An International Common Mechanism for DNA Synthesis Screening"),
    keywords = "",
    url = "https://github.com/DNA-Screening/CommonMechanism",
    packages = setuptools.find_packages(where="common_mechanism", exclude=("tests",)),
    package_data = { 'src/databases': ['*'] },
    install_requires=["Biopython", "pandas", "plotly", "matplotlib", "kaleido", "pybedtools", "pytaxonkit"],
    long_description=read('README.md'),
)
