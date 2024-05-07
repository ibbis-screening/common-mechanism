# commec: an open-source, globally available tool for DNA sequence screening

The `commec` package is a free and open-source tool for DNA sequence screening that is part of the
[Common Mechanism for DNA Synthesis screening](https://ibbis.bio/common-mechanism/).

Introduction
============
The Common Mechanism offers three sub-commands through the `commec` entrypoint:

    screen  Run Common Mechanism screening on an input FASTA.
    flag    Parse all .screen files in a directory and create two CSVs file of flags raised
    split   Split a multi-record FASTA file into individual files, one for each record

The `screen` command runs an input FASTA through four steps:

  1. Biorisk scan (uses a hmmer search against custom databases)
  2. Regulated protein scan (uses a BLASTX or DIAMOND search against NCBI nr)
  3. Regulated nucleotide scan (uses BLASTN against NCBI nt)
  4. Benign scan (users hmmer, cmscan and BLASTN against custom databases)

The `.screen` file produced by that pipeline can be passed to `flag` to produce two output CSVs.
`flags.csv` will have the following columns:

    filename:                 .screen file basename
    biorisk:                  "F" if flagged, "P" if no flags
    virulence_factor:         "F" if flagged, "P" if no flags
    regulated_virus:          "F" if flagged, "P" if no flags, "Err" if error logged
    regulated_bacteria:       "F" if flagged, "P" if no flags, "Err" if error logged
    regulated_eukaryote:      "F" if flagged, "P" if no flags, "Err" if error logged
    mixed_regulated_non_reg:  "F" if flagged, "P" if no flags, "Err" if error logged
    benign:                   "F" if not cleared, "P" if all cleared, "-" if not run

The flags_recommended CSV just has two columns, "filename" and "recommend_flag_or_pass". The
recommendation is based on the following decision flow:



Documentation
=============
The online documentation is located at the
[GitHub Wiki](https://github.com/ibbis-screening/common-mechanism/wiki).

Development
=======
Development dependencies are managed through a conda environment. Install conda, then make sure
that [your channels are configured correctly](http://bioconda.github.io/).

```
conda env create -f environment.yml
conda activate commec-dev
```

From here, you should have an interactive version of the package installed via (`pip -e .`) as well
as the necessary shell dependencies.

About
=====
The Common Mechanism is a project of [IBBIS](https://ibbis.bio), the International Biosecurity and
Biosafety Initiative for Science. From 2021-2023, the software and databases were developed by a
team of technical consultants working with the Nuclear Threat Initiative, led by Dr. Nicole Wheeler
of the University of Birmingham, and including contributions from Brittany Rife Magalis of the
University of Louisville and Jennifer Lu of the Center for Computational Biology at Johns Hopkins
University. In 2024, IBBIS became the home of the project
