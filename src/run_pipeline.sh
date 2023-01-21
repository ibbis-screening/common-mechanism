#! /usr/bin/env bash

##############################################################################
#run_pipeline.sh runs the Common Mechanism against a specified QUERY file.
#
#Copyright (C) 2022-2023 NTI|Bio 
#This file is part of the CommonMechanism 
##############################################################################
# Usage:
#   src/run_pipeline.sh
#       -q INPUT_FILE
#       -d DATBASE_FOLDER 
#       -t THREADS (default: 1) 
#       -o OUTPUT (output prefix, default: query name) 
#   Optional parameters:
#       -c = clean up intermediate files (default: no cleanup)
#       -b = use blast instead of diamond (default: diamond) 
# Example Usage: src/run_pipeline.sh-q test_folder/[$name].fasta -d databases/ -p 5 -t 1 -o out_prefix 
##############################################################################
# set parameters for the run
#set -eu 
PROCESSES=5         #number of processes to run at once
THREADS=1           #threads per process
QUERY=""            #query input file 
OUTPUT=""           #output prefix
CLEANUP=0           #cleanup files (default: false)" 
BLAST=0             #blast or diamond (default: diamond) 
DB_PATH="" 

function print_usage() {
    echo ""
    echo " Usage: src/run_pipeline.sh -q QUERY -d DB_PATH/ -o OUTPUT [-t THREADS]"
    echo "    QUERY           query file to align to each database (required)"
    echo "    OUTPUT          output prefix for alignments (default: query prefix)"
    echo "    DB_PATH         path to detabases (required)"
    echo "    THREADS         number of threads for each database run (default: 1)"
    echo " OPTIONAL FLAGS" 
    echo "    -c              tidy up intermediate screening files afterward"
    echo "    -b              run blast for protein screen [default: diamond]"
    
}

#Get options from user
while getopts "t:d:q:o:c:b:" OPTION
    do
        case $OPTION in
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
                CLEANUP=1
                ;;
            b)
                BLAST=1
                ;;
            \?)
                echo "Usage: src/run_pipeline.sh -q QUERY -d DB_PATH/ -o OUTPUT [-t THREADS]"
                echo "  QUERY           query file to align to each database (required)"
                echo "  OUTPUT          output prefix for alignments (default: query prefix)"
                echo "  DB_PATH         path to detabases (required)"
                echo "  THREADS         number of threads for each database run (default: 1)"
                echo "OPTIONAL FLAGS" 
                echo "  -c              tidy up intermediate screening files afterward"
                echo "  -b              run blast for protein screen [default: diamond]"
                exit
                ;;
        esac
    done

#INITIALIZING: Check for values
#Check input query 
echo " >> STEP 0: Checking for valid options..." 
if [ "$QUERY" == "" ]; then
    echo " ERROR: no input query file specified"
    print_usage
    exit 1 
elif [ ! -f "$QUERY" ]; then 
    echo " ERROR: specified input query file does not exist"
    print_usage
    exit 1
fi
#If output not specified, set value 
if [ "$OUTPUT" == "" ]; then
    OUTPUT=${QUERY//.fasta/}
    echo "Output handle = $OUTPUT"
fi
#Check input database folder 
if [ "$DB_PATH" == "" ]; then
    echo " ERROR: screening database path not specified"
    print_usage
    exit 1
elif [ ! -d "$DB_PATH" ]; then
    echo " ERROR: specified database folder does not exist"
    print_usage
    exit 1
else
    #Test database downloads 
    biorisk_test=${DB_PATH}/biorisk_db/biorisk.hmm
    benign_test=${DB_PATH}/benign_db/benign.hmm
    split_nr_test=${DB_PATH}/nr_dmnd/nr.1.dmnd
    if [ ! -f "$biorisk_test" ]; then
        echo " ERROR: database folder does not contain biorisk database"
        exit 1 
    fi 
    if [ ! -f "$benign_test" ]; then
        echo " ERROR: database folder does not contain benign database" 
        exit 1 
    fi 
    if [ ! -f "$split_nr_test" ]; then 
        echo " ERROR: database folder does not contain diamond database"
        exit 1
    fi 
fi
# export DB_PATH=${DB_PATH}

#################################################################
CM_DIR="$( dirname "$0" )"

start_time=$(date)
echo -e " >> STARTED AT $start_time"

# Step 1: biorisk DB scan
echo " >> STEP 1: Running biorisk hmm scan..." 
echo -e "\t...running transeq" 
transeq $QUERY ${OUTPUT}.faa -frame 6 -clean &>> ${OUTPUT}.tmp
echo -e "\t...running hmmscan" 
hmmscan --domtblout ${OUTPUT}.biorisk.hmmsearch ${DB_PATH}/biorisk_db/biorisk.hmm ${OUTPUT}.faa &>> ${OUTPUT}.tmp
echo -e "\t...checking hmmscan results"
python ${CM_DIR}/check_biorisk.py -i ${OUTPUT}.biorisk.hmmsearch --database ${DB_PATH}/biorisk_db/

s1_time=$(date)
echo -e "    STEP 1 completed at $s1_time\n"

# Step 2: taxon ID/protein screening
echo " >> STEP 2: Checking regulated pathogen proteins..."
if [ "$BLAST" = 1 ]; then
    if [ ! -f "${OUTPUT}".nr.blastx ]; then
        echo -e "\t...running run_blastx.sh"
        ${CM_DIR}/run_blastx.sh -d $DB_PATH/nr_blast/nr -q $QUERY -o ${OUTPUT}.nr -t $THREADS
    fi
    echo -e "\t...checking blast results"
    python ${CM_DIR}/check_reg_path.py -i ${OUTPUT}.nr.blastx --benign-db $DB_PATH/benign_db/ --biorisk-db $DB_PATH/biorisk_db/
else 
    if [ ! -f "${OUTPUT}".nr.dmnd ]; then
       echo -e "\t...running run_diamond.sh"
        ${CM_DIR}/run_diamond.sh -d $DB_PATH/nr_dmnd/ -i $QUERY -o ${OUTPUT}.nr -t $THREADS -p $PROCESSES 
        # cat ${OUTPUT}.nr* > ${OUTPUT}.nr.dmnd
    fi
    echo -e "\t...checking diamond results"
    python ${CM_DIR}/check_reg_path.py -i ${OUTPUT}.nr.dmnd --benign-db $DB_PATH/benign_db/ --biorisk-db $DB_PATH/biorisk_db/
fi

s2_time=$(date)
echo -e "    STEP 2 completed at $s2_time\n"

# nucleotide screening
echo " >> STEP 3: Checking regulated pathogen nucleotides..."
echo -e "\t...fetching noncoding regions"
if [ "$BLAST" = 1 ]; then
    python ${CM_DIR}/fetch_nc_bits.py ${OUTPUT}.nr.blastx ${QUERY}
else
    python ${CM_DIR}/fetch_nc_bits.py ${OUTPUT}.nr.dmnd ${QUERY}
fi

if [ -f "${OUTPUT}".noncoding.fasta ]; then 
    echo -e "\t...running blastn"
    blastn -query ${OUTPUT}.noncoding.fasta -db ${DB_PATH}/nt_blast/nt -out ${OUTPUT}.nt.blastn -outfmt "7 qacc stitle sacc staxids evalue bitscore pident qlen qstart qend slen sstart send" -max_target_seqs 50 -num_threads 8 -culling_limit 5 -evalue 10
    echo -e "\t...checking blastn results"
    python ${CM_DIR}/check_reg_path.py -i ${OUTPUT}.nt.blastn --benign-db $DB_PATH/benign_db/ --biorisk-db $DB_PATH/biorisk_db/
else 
    echo "\t...skipping nucleotide search"
fi

s3_time=$(date)
echo -e "    STEP 3 completed at $s3_time\n"

# Step 3: benign DB scan
#date
echo -e " >> STEP 4: Checking any pathogen regions for benign components..."
hmmscan --domtblout ${OUTPUT}.benign.hmmsearch ${DB_PATH}/benign_db/benign.hmm ${OUTPUT}.faa &>>${OUTPUT}.tmp
blastn -db ${DB_PATH}/benign_db/benign.fasta -query $QUERY -out ${OUTPUT}.benign.blastn -outfmt "7 qacc stitle sacc staxids evalue bitscore pident qlen qstart qend slen sstart send" -evalue 1e-5

python ${CM_DIR}/check_benign.py -i ${OUTPUT} --sequence ${QUERY} -d ${DB_PATH}/benign_db/ # added the original file path here to fetch sequence length, can tidy this

finish_time=$(date)
echo -e " >> COMPLETED AT $finish_time"

# Visualising outputs; functional characterization

#python ${CM_DIR}/viz_outputs.py ${OUTPUT} # turning off until file write permissions are resolved

if [ "$CLEANUP" == 1 ]
then
    if [ -f "${OUTPUT}".reg_path_coords.csv ]
    then
        rm ${OUTPUT}.reg_path_coords.csv
    fi
    rm ${OUTPUT}.*hmmsearch ${OUTPUT}.*blastx ${OUTPUT}.*blastn
fi

rm ${OUTPUT}*.tmp
