#! /usr/bin/env bash

# usage: src/run_pipeline.sh test_folder/[$name].fasta

#PKG_HOME=`which run_pipeline.sh | sed 's/\/src\/run_pipeline.sh//g'`

# run a config file that sets relevant paths
dirname = "$( dirname "$0" )"
source ${dirname}/config

query=$1 # the file name
name=${query//*\//} # strip out any directory info
name=${name//.fasta/} # the prefix detailing the name of the sequence

# Step 1: biorisk DB scan
transeq $query ${name}.faa -frame 6 -clean &>/dev/null
hmmscan --domtblout ${name}.biorisk.hmmsearch biorisk/biorisk.hmm ${name}.faa &>/dev/null
python -m check_biorisk ${name}

# Step 2: taxon ID
#if ! [ -e "${name}.nr.blastx" ]; # if the file already exists, don't remake it - can remove this upon release (if added back, remember then before blastx command)
#blastx -db nr -query $query -out ${name}.nr.blastx -outfmt "7 qacc stitle sacc staxids evalue bitscore pident qlen qstart qend slen sstart send" -max_target_seqs 500 -num_threads 8 -culling_limit 5 -evalue 1e-10 -word_size 6 -threshold 21 -window_size 40 -matrix BLOSUM62 -gapopen 11 -gapextend 1 -seg yes
#fi

src/run_blastx.sh -d nr -q $query -o ${name}.nr

### IF A HIT TO A REGULATED PATHOGEN, PROCEED, OTHERWISE CAN FINISH HERE ONCE TESTING IS COMPLETE ####
python -m check_reg_path ${name}

# Step 3: benign DB scan
hmmscan --domtblout ${name}.benign.hmmsearch benign/benign.hmm ${name}.faa &>/dev/null
blastn -db benign/benign.fasta -query $query -out ${name}.benign.blastn -outfmt "7 qacc stitle sacc staxids evalue bitscore pident qlen qstart qend slen sstart send" -evalue 1e-5

python -m check_benign ${name} ${query} # added the original file path here to fetch sequence length, can tidy this

# Visualising outputs; functional characterization

#python -m viz_outputs ${name} # turning off until file write permissions are resolved

#rm ${query}.reg_path_coords.csv $name.*hmmsearch $name.*blastx $name.*blastn
