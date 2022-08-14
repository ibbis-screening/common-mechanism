#! /usr/bin/env bash

# usage: src/run_pipeline.sh test_folder/[$name].fasta

# need to make this flexible re: where people store their databases
export PYTHONPATH=$PYTHONPATH/src
export PFAMDB=databases
export BLASTDB=$BLASTDB:./databases

query=$1 # the file name
name=${query//*\//} # strip out any directory info
name=${name//.fasta/} # the prefix detailing the name of the sequence

# Step 1: biorisk DB scan
transeq $query ${name}.faa -frame 6 -clean &>/dev/null
hmmscan --domtblout ${name}.biorisk.hmmsearch biorisk/biorisk.hmm ${name}.faa &>/dev/null
python src/check_biorisk.py ${name}

# Step 2: taxon ID
if ! [ -e "${name}.nr.diamond" ];
then src/run_diamond.sh -d ../databases/nr_dmnd/ -i ${query}.fasta -o ${name}.nr.diamond -t 4 -p 4
cat $name.*.tsv > $name.nr.diamond
rm $name.*.tsv
fi

# here need to flag any uncovered regions >200 bp and blastn them
# currently running this as a step in check_reg_path

### IF A HIT TO A REGULATED PATHOGEN, PROCEED, OTHERWISE CAN FINISH HERE ONCE TESTING IS COMPLETE ####
python src/check_reg_path_diamond_proteins.py ${name}

# nucleotide screening

python src/fetch_nc_bits.py ${name} ${query}
blastn -query ${name}_nr.fasta -db nt

# Step 3: benign DB scan
hmmscan --domtblout ${name}.benign.hmmsearch benign/benign.hmm ${name}.faa &>/dev/null
blastn -db benign/benign.fasta -query $query -out ${name}.benign.blastn -outfmt "7 qacc stitle sacc staxids evalue bitscore pident qlen qstart qend slen sstart send" -evalue 1e-5

python src/check_benign.py ${name} ${query} # added the original file path here to fetch sequence length, can tidy this

# Visualising outputs; functional characterization

if [ ! -d //figures ]; then
  mkdir -p figures;
fi
python src/viz_outputs.py ${name}

#rm ${query}.reg_path_coords.csv $name.*hmmsearch $name.*blastx $name.*blastn
