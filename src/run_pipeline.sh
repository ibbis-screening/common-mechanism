#! /usr/bin/env bash

##############################################################################
#run_pipeline.sh runs the Common Mechanism against a specified QUERY file.
##############################################################################
# Usage: src/run_pipeline.sh test_folder/[$name].fasta

# parameters for the run
set -eu
PROCESSES=5         #number of processes to run at once
THREADS=1           #threads per process
DB_PATH=""
QUERY=""
OUTPUT=""
CLEANUP=1

#Get options from user
while getopts "p:t:d:q:o:c:" OPTION
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
            q)
                QUERY=$OPTARG
                ;;
            o)
                OUTPUT=$OPTARG
                ;;
            c)
                CLEANUP=$OPTARG
                ;;
            \?)
                echo "Usage: src/run_pipeline.sh -d DB_PATH -q QUERY -s OUTPUT [-p PROCESSES -t THREADS]"
                echo "  DB_PATH              location (folder) of database (required)"
                echo "  QUERY           query file to align to each database (required)"
                echo "  OUTPUT          output prefix for alignments (default: query prefix)"
                echo "  PROCESSES       number of databases to evaluate (default: 5)"
                echo "  THREADS         number of threads for each database run (default: 1)"
                echo "  CLEANUP         tidy up intermediate screening files afterward?"
                exit
                ;;
        esac
    done

#Check for values
if [ "$DB_PATH" == "" ] && [ "$QUERY" == "" ]
then
    echo "Usage: run_blastx.sh -d DB_PATH -q QUERY -s OUTPUT [-p PROCESSES -t THREADS]"
        echo "  DB_PATH         location (folder) of database (required)"
        echo "  QUERY           query file to align to each database (required)"
        echo "  OUTPUT          output prefix for alignments (default: out)"
        echo "  PROCESSES       number of databases to evaluate (default: 5)"
        echo "  THREADS         number of threads for each database run (default: 1)"
        echo "  CLEANUP         tidy up intermediate screening files afterward? 0 = no, 1 = y"
    exit
fi


if [ "$OUTPUT" == "" ]
then
    OUTPUT=$(basename "$QUERY" | cut -d. -f1)
fi

#Check for database
echo " >> Checking for valid options..."
if [ -d "$DB_PATH" ]
then
    #Directory exists, check for at least one blastx db file
    if [ ! -f "$DB_PATH"/*.pal ]
    then
        echo " ERROR: no blast database in $DB_PATH"
        exit
    fi
else
    echo " ERROR: database folder $DB_PATH does not exist"
    exit
fi

#Check for input file
if [ ! -f  $QUERY ]
then
    echo " ERROR: input file $QUERY does not exist"
    exit
fi

#################################################################

# run a config file that sets relevant paths
dirname="$( dirname "$0" )"

export CM_DIR=$dirname
export PYTHONPATH=$dirname

source ${CM_DIR}/../config

name=${QUERY//*\//} # strip out any directory info
name=${QUERY//.fasta/} # the prefix detailing the name of the sequence

# Step 1: biorisk DB scan
echo " >> Running biorisk HMM scan..."
transeq $QUERY ${name}.faa -frame 6 -clean &>tmp
hmmscan --domtblout ${name}.biorisk.hmmsearch biorisk/biorisk.hmm ${name}.faa &>tmp
python ${CM_DIR}/check_biorisk.py ${name}

# Step 2: taxon ID
# protein screening
echo " >> Running taxid screen for regulated pathogens..."
${CM_DIR}/run_blastx.sh -d ${DB_PATH}/nr -q $QUERY -o ${name}.nr -t $THREADS

python ${CM_DIR}/check_reg_path_proteins.py ${name}

# nucleotide screening

#python src/fetch_nc_bits.py ${name} ${QUERY}
#blastn -query ${name}_nc.fasta -db ${DB_PATH}/../nt_blast/nt

# .nt.blastn -outfmt "7 qacc stitle sacc staxids evalue bitscore pident qlen qstart qend slen sstart send" -max_target_seqs 500 -num_threads 8 -culling_limit 5 -evalue 30 -word_size 7 -window_size 40 -gapopen 11 -gapextend 1'

#python ${CM_DIR}/check_reg_path_nt.py ${name}

### IF A HIT TO A REGULATED PATHOGEN, PROCEED, OTHERWISE CAN FINISH HERE ONCE TESTING IS COMPLETE ####

# Step 3: benign DB scan
echo " >> Checking any pathogen regions for benign components..."
hmmscan --domtblout ${name}.benign.hmmsearch benign/benign.hmm ${name}.faa &>/dev/null
blastn -db ${PFAMDB}/benign/benign.fasta -query $QUERY -out ${name}.benign.blastn -outfmt "7 qacc stitle sacc staxids evalue bitscore pident qlen qstart qend slen sstart send" -evalue 1e-5

python ${CM_DIR}/check_benign.py ${name} ${QUERY} # added the original file path here to fetch sequence length, can tidy this

echo " >> Done with screening!"

# Visualising outputs; functional characterization

#python ${CM_DIR}/viz_outputs.py ${name} # turning off until file write permissions are resolved

if [ "$CLEANUP" == 1 ]
then
    if [ -f "${name}".reg_path_coords.csv ]
    then
        rm ${name}.reg_path_coords.csv
    fi
    rm ${name}.*hmmsearch ${name}.*blastx ${name}.*blastn
fi

rm tmp
