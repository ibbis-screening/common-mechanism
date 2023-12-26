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
    license = "MIT",
    keywords = "",
    url = "https://github.com/nwheeler443/CommonMechanism",
    packages = setuptools.find_packages(where="common_mechanism", exclude=("tests",)),
    package_data = { 'src/databases': ['*'] },
    install_requires=["Biopython", "pandas", "taxoniq", "gdown", "plotly", "matplotlib", "kaleido", "pybedtools"],
    long_description=read('README.md'),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Topic :: Scientific/Engineering :: Bio-Informatics"
        "License :: OSI Approved :: MIT License",
    ],
)
