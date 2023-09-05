# An International Common Mechanism for DNA Synthesis Screening
CommonMechanism is designed for screening DNA for regulated pathogens that 
pose a significant biorisk. In order to run this pipeline, please first see [Installation](#installation) and [Required Data](#required-data) for setup. 

Table of Contents:
1. [Installation](#installation) 
2. [Required Data](#required-data) 
3. [Running the CommonMechanism](#running-the-commonmechanism)
4. [User Survey](#user-survey)
5. [Author Information](#author-information)

# Installation 
In addition to the CommonMechanism source code, users may need to install a number of software packages that will be called by the CommonMechanism. This section details the various dependencies

## HMMER search
The CommonMechanism pipeline requires:
 * hmmer 
 * transeq (from the EMBOSS package)

to run the biorisk search.

## Python Packages 
The CommonMechanism pipeline has the following python package dependencies:
 * pandas
 * plotly 
 * taxoniq 
 * ncbi-taxon-db
 * biopython
 * matplotlib
 * kaleido
 * bedtools
 * pybedtools 

## BLAST Aligner 
The CommonMechanism requires the BLAST code executable which is found at https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/. 

For more details on BLAST, see https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html. 

## DIAMOND Aligner 
The CommonMechanism pipeline requires users to install [DIAMOND](https://github.com/bbuchfink/diamond "DIAMOND github") by the following commands: 
   
    wget https://github.com/bbuchfink/diamond/archive/v2.0.13.tar.gz
    tar xzf v2.0.13.tar.gz
    mkdir diamond-2.0.13/bin
    cd diamond-2.0.13/bin
    make -j8
    sudo make install

Following installation of all required packages, please see [Required Data](#required-data)

# RNA family search

The CommonMechanism uses Infernal to search sequences for known benign RNA families. It can be installed using:

`conda install -c bioconda infernal`

# Required Data 
The following databases will need to be installed prior to running the CommonMechanism pipeline. Please place all databases in the same folder. The resulting file structure (and required storage space) will be as follows:

    databases/
    databases/nr_dmnd (~276 Gb) 
    databases/nr_blast (~375 Gb)   (optional)
    databases/nt_blast (~280 Gb) 
    databases/benign_db (~2.2 Gb)
    databases/biorisk_db (~1.0 Gb)

## Benign and Biorisk Databases 
To download the benign_db and biorisk_db files, please contact Nicole Wheeler (n.wheeler@bham.ac.uk) and she will send you a download link.

## BLAST Database
The Common Mechanism requires the BLAST nt database. The files for the BLAST nt database are located at https://ftp.ncbi.nlm.nih.gov/blast/db/. Users may download the BLAST database (including 700+ files) within the `nt_blast/` folder by using the following command:

      update_blastdb.pl --passive --decompress nt
   
This command requires the BLAST executables to be installed (see [BLAST Download](#blast-aligner))
Users can also download the BLAST nr database to check the relative speed of BLAST vs DIAMOND search on their machine:

      update_blastdb.pl --passive --decompress nr

## DIAMOND Database
The CommonMechanism database will be provided via Amazon AWS. Please email Jennifer Lu (jennifer.lu717@gmail.com) for instructions to download the database. 

# Running the Common Mechanism 
The basic command line for running the Common Mechanism is as follows:

      /run_pipeline.sh -q ${QUERY} -o ${OUTPUT} -d ${DATABASE FOLDER}

The following required parameters must be specified:
- `${QUERY}` - The sequence file to test
- `${OUTPUT}` - User-specified output prefix (All output files will begin with ${OUTPUT})
- `${DATABASE FOLDER}` - The path to the `databases/` folder setup as described in [Required Data](#required-data)

# User survey
https://docs.google.com/forms/d/1LqzEH3XkFHBtMHNzGW9i4vRU7TyabBnNIDOFQBXa-4s/edit?usp=sharing

# Author information 
Nicole Wheeler (n.wheeler@bham.ac.uk)
Jennifer Lu (jennifer.lu717@gmail.com)

Last Updated On: 11/10/2022
