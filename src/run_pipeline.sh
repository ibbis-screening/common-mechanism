#!/bin/bash

<<<<<<< HEAD
export PFAMDB=./databases/
export BLASTDB=$BLASTDB:../
=======
# usage: src/run_pipeline.sh test_sequences/[$name].fasta

export PFAMDB=$PFAMDB:./databases
export BLASTDB=$BLASTDB:./databases
>>>>>>> cfc60b3ce8102364e0710d1e86026fbb25e77671

query=$1 # the file name
name=${query//*\//}
name=${name//.fasta/} # the prefix detailing the name of the sequence

# biorisk DB scan
transeq $query ${name}.faa -frame 6 -clean &>/dev/null
hmmscan --domtblout ${name}.biorisk.hmmsearch biorisk/biorisk.hmm ${name}.faa &>/dev/null
python src/check_biorisk.py ${name}

# taxon ID
# blastx -db nr -query $query -out ${name}.nr.blastx -outfmt "7 qacc stitle sacc staxids evalue bitscore pident qlen qstart qend slen sstart send" -max_target_seqs 500 -evalue 1e-5 -remote
# not enough room on computer
blastx -db nr -query $query -out ${name}.nr.blastx -outfmt "7 qacc stitle sacc staxids evalue bitscore pident qlen qstart qend slen sstart send" -max_target_seqs 500 -evalue 1e-5 -num_threads 8

### IF A HIT TO A REGULATED PATHOGEN, PROCEED, OTHERWISE CAN FINISH HERE ONCE TESTING IS COMPLETE ####
python src/check_reg_path.py ${name}

# benign DB scan
hmmscan --domtblout ${name}.benign.hmmsearch benign/benign.hmm ${name}.faa &>/dev/null
<<<<<<< HEAD
blastn -db benign.fasta -query $query -out ${name}.benign.blastn -outfmt "7 qacc stitle sacc staxids evalue bitscore pident qlen qstart qend slen sstart send" -evalue 1e-5
=======
blastn -db benign/benign.fasta -query $query -out ${name}.benign.blastn -outfmt "7 qacc stitle sacc staxids evalue bitscore pident qlen qstart qend slen sstart send" -evalue 1e-5
>>>>>>> cfc60b3ce8102364e0710d1e86026fbb25e77671

python src/check_benign.py ${name} ${query} # added the original file path here to fetch sequence length, can tidy this

# functional characterization

<<<<<<< HEAD
# rm ${name}.reg_path_coords.csv $name.*hmmsearch $name.*blastx $name.*blastn
=======
python src/viz_outputs.py ${name}

#rm ${name}.reg_path_coords.csv $name.*hmmsearch $name.*blastx $name.*blastn
>>>>>>> cfc60b3ce8102364e0710d1e86026fbb25e77671
