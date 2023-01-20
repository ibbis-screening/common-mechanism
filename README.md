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
The CommonMechanism pipeline has the following package dependencies:
 * plotly 
 * taxoniq 
 * biopython
 * matplotlib
 * kaleido
 * pybedtools 

Additionally, the CommonMechanism pipeline requires users to install [DIAMOND](https://github.com/bbuchfink/diamond "DIAMOND github") by the following commands: 
   
    wget https://github.com/bbuchfink/diamond/archive/v2.0.13.tar.gz
   tar xzf v2.0.13.tar.gz
   mkdir diamond-2.0.13/bin
   cd diamond-2.0.13/bin
   cmake -DWITH_ZSTD=ON ZSTD_LIBRARY=~/miniconda3/lib -DBLAST_INCLUDE_DIR=$HOME/BLAST2.11/include/ncbi-tools++ -DBLAST_LIBRARY_DIR=$HOME/BLAST2.11/lib ..
   make -j8
   sudo make install

Following installation of all required packages, please see [Required Data](#required-data)

# Required Data 
The following databases will need to be installed prior to running the CommonMechanism pipeline. Please place all databases in the same folder. The resulting file structure will be as follows:

    databases/
    databases/split_nr
    databases/benign_db
    databases/biorisk_db 

## Benign and Biorisk Databases 
https://drive.google.com/uc?id=1INLOrC_vLMxEmtFCW4rzyZe5z6hWD1Qu

# Running the Common Mechanism 


# User survey
https://docs.google.com/forms/d/1LqzEH3XkFHBtMHNzGW9i4vRU7TyabBnNIDOFQBXa-4s/edit?usp=sharing

# Author information 
Nicole Wheeler (n.wheeler@bham.ac.uk)
Jennifer Lu (jennifer.lu717@gmail.com)

Last Updated On: 11/10/2022
