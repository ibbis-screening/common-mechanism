#! /usr/bin/env bash

#####################################################################
#run_diamond.sh runs DIAMOND against a specified database.
#   Buchfink B, Reuter K, Drost HG, "Sensitive protein alignments at 
#   tree-of-life scale using DIAMOND", 
#   Nature Methods 18, 366–368 (2021). 
#   doi:10.1038/s41592-021-01101-x
#Script Author: Jennifer Lu, jennifer.lu717@gmail.com
#Updated: 05/25/2022 
##################################################################### 
#Usage: run_diamond.sh -d MY_DB -i INPUT_FILE -o OUTPUT_FILE [-p PROCESSES -t THREADS] 

#set -eux #debug mode
set -eu
PROCESSES=5         #number of processes to run at once 
THREADS=1           #threads per process  
DB_PATH=""
INPUT=""
OUTPUT="out" 

#Get options from user
while getopts "p:t:d:i:o:" OPTION
    do 
        case $OPTION in
            p) 
                PROCESSES=$OPTARG
                ;;
            t) 
                THREADS=$OPTARG
                ;;
            d) 
                DB_PATH=$OPTARG
                ;;
            i) 
                INPUT=$OPTARG
                ;;
            o) 
                OUTPUT=$OPTARG
                ;;
            \?)
                echo "Usage: run_diamond.sh -d MY_DB -i INPUT_FILE -o OUTPUT_FILE [-p PROCESSES -t THREADS]"
                echo "  MY_DB           location of database (required)"
                echo "  INPUT_FILE      input file to align to each database (required)" 
                echo "  OUTPUT_FILE     output prefix for alignments (default: out)" 
                echo "  PROCESSES       number of databases to evaluate (default: 5)"
                echo "  THREADS         number of threads for each database run (default: 1)"
                exit 
                ;;  
        esac 
    done

#Check for values
if [ "$DB_PATH" == "" ] && [ "$INPUT" == "" ]
then
    echo "Usage: run_diamond.sh -d MY_DB -i INPUT_FILE -o OUTPUT_FILE [-p PROCESSES -t THREADS]"
    echo "  MY_DB           location of database (required)"
    echo "  INPUT_FILE      input file to align to each database (required)" 
    echo "  OUTPUT_FILE     output prefix for alignments (default: out)" 
    echo "  PROCESSES       number of databases to evaluate (default: 5)"
    echo "  THREADS         number of threads for each database run (default: 1)"
    exit 
fi
#Check for database 
echo " >> Checking for Valid Options..." 
if [ -d $DB_PATH ]
then 
    #Directory exists, check for at least one diamond db file 
    if [ ! -f $DB_PATH/*.1.dmnd ]
    then 
        echo " ERROR: diamond database $DB_PATH/nr.1.dmnd does not exist"
        exit
    fi
else
    echo " ERROR: nr diamond database folder $DB_PATH does not exist" 
    exit
fi
  
#Check for input file 
if [ ! -f  $INPUT ] 
then
    echo " ERROR: input file $INPUT does not exist" 
    exit
fi      
  
ls ${DB_PATH}/*.dmnd | parallel --will-cite -j ${PROCESSES} diamond blastx -d ${DB_PATH}/uniref90.{%}.dmnd --fast --threads ${THREADS} -q ${INPUT} -o ${OUTPUT}.{%}.tsv --frameshift 15 --range-culling --outfmt 6 qseqid stitle sseqid staxids evalue bitscore pident qlen qstart qend slen sstart send

