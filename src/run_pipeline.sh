#! /usr/bin/env bash

##############################################################################
#run_pipeline.sh runs the Common Mechanism against a specified QUERY file.
##############################################################################
# Usage: src/run_pipeline.sh test_folder/[$name].fasta

# parameters for the run
set -eu
PROCESSES=5         #number of processes to run at once
THREADS=1           #threads per process
QUERY=""
OUTPUT=""
CLEANUP=0
DB_PATH=""

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
                echo "Usage: src/run_pipeline.sh -q QUERY -o OUTPUT [-p PROCESSES -t THREADS]"
                echo "  QUERY           query file to align to each database (required)"
                echo "  OUTPUT          output prefix for alignments (default: query prefix)"
                echo "  DB_PATH         path to detabases (required)"
                echo "  PROCESSES       number of databases to evaluate (default: 5)"
                echo "  THREADS         number of threads for each database run (default: 1)"
                echo "  CLEANUP         tidy up intermediate screening files afterward?"
                exit
                ;;
        esac
    done

#Check for values
if [ "$QUERY" == "" ]
then
    echo "Usage: src/run_pipeline.sh -q QUERY -o OUTPUT [-p PROCESSES -t THREADS]"
        echo "  QUERY           query file to align to each database (required)"
        echo "  OUTPUT          output prefix for alignments (default: query prefix)"
        echo "  DB_PATH         path to detabases (required)"
        echo "  PROCESSES       number of databases to evaluate (default: 5)"
        echo "  THREADS         number of threads for each database run (default: 1)"
        echo "  CLEANUP         tidy up intermediate screening files afterward? 0 = no, 1 = y"
    exit
fi


if [ "$OUTPUT" == "" ]
then
    OUTPUT=$(basename "$QUERY" | cut -d. -f1)
fi

echo " >> Checking for valid options..."

if [ "$DB_PATH" == "" ]
then
    echo "Please specify the path to screening databases"
fi
export $DB_PATH

#Check for input file
if [ ! -f  $QUERY ]
then
    echo " ERROR: input file $QUERY does not exist"
    exit
fi

#################################################################

dirname="$( dirname "$0" )"

export CM_DIR=$dirname
export PYTHONPATH=$dirname

# Step 1: biorisk DB scan
date
echo " >> Running biorisk HMM scan..."
transeq $QUERY ${OUTPUT}.faa -frame 6 -clean &>> ${OUTPUT}.tmp
hmmscan --domtblout ${OUTPUT}.biorisk.hmmsearch -E 1e-10 ${DB_PATH}/biorisk/biorisk.hmm ${OUTPUT}.faa &>> ${OUTPUT}.tmp
python ${CM_DIR}/check_biorisk.py ${OUTPUT}

# Step 2: taxon ID
# protein screening
date
echo " >> Running taxid screen for regulated pathogen proteins..."
${CM_DIR}/run_blastx.sh -d $DB_PATH/nr_blast/nr -q $QUERY -o ${OUTPUT}.nr -t $THREADS

python ${CM_DIR}/check_reg_path_prot.py ${OUTPUT}

# nucleotide screening
date
echo " >> Running taxid screen for regulated pathogen nucleotides..."
python ${CM_DIR}/fetch_nc_bits.py ${OUTPUT} ${QUERY}
if [ -f "${OUTPUT}"_nc.fasta ]
then blastn -query ${OUTPUT}_nc.fasta -db ${DB_PATH}/nt_blast/nt -out ${OUTPUT}.nt.blastn -outfmt "7 qacc stitle sacc staxids evalue bitscore pident qlen qstart qend slen sstart send" -max_target_seqs 50 -num_threads 8 -culling_limit 5 -evalue 10
python ${CM_DIR}/check_reg_path_nt.py ${OUTPUT}
fi

# Step 3: benign DB scan
date
echo " >> Checking any pathogen regions for benign components..."
hmmscan --domtblout ${OUTPUT}.benign.hmmsearch -E 1e-10 ${DB_PATH}/benign/benign.hmm ${OUTPUT}.faa &>>${OUTPUT}_tmp
blastn -db ${DB_PATH}/benign/benign.fasta -query $QUERY -out ${OUTPUT}.benign.blastn -outfmt "7 qacc stitle sacc staxids evalue bitscore pident qlen qstart qend slen sstart send" -evalue 1e-5

python ${CM_DIR}/check_benign.py ${OUTPUT} ${QUERY} # added the original file path here to fetch sequence length, can tidy this

date
echo " >> Done with screening!"

# Visualising outputs; functional characterization

#python ${CM_DIR}/viz_outputs.py ${OUTPUT} # turning off until file write permissions are resolved

if [ "$CLEANUP" == 1 ]
then
    if [ -f "${OUTPUT}".reg_path_coords.csv ]
    then
        rm ${OUTPUT}.reg_path_coords.csv
    fi
    if [ -f "${OUTPUT}".reg_path_coords_nt.csv ]
    then
        rm ${OUTPUT}.reg_path_coords_nt.csv
    fi
    rm ${OUTPUT}.*hmmsearch ${OUTPUT}.*blastx ${OUTPUT}.*blastn
fi

rm ${OUTPUT}*.tmp
