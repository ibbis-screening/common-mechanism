#!/bin/bash

query=$1 # the file prefix before. fasta
echo $query

blastn -db nt -query ../queries/${query}.fasta -out blast/${query}.nr.blastn -outfmt "7 qacc stitle sacc staxids evalue bitscore pident qlen qstart qend slen sstart send" -max_target_seqs 500 -evalue 10

blastx -db swissprot -query ../queries/${query}.fasta -out blast/${query}.sprot.blastx -outfmt "7 qacc stitle sacc staxids evalue bitscore pident qlen qstart qend slen sstart send" -max_target_seqs 500 -evalue 1e-5

blastx -db VFDB_setA_pro.fas -query ../queries/${query}.fasta -out blast/${query}.vfdb.blastx -outfmt "7 qacc stitle sacc staxids evalue bitscore pident qlen qstart qend slen sstart send" -max_target_seqs 1000 -evalue 10

blastn -db blacklist.fasta -query ../queries/${query}.fasta -out blast/${query}.blacklist.blastn -outfmt "7 qacc stitle sacc staxids evalue bitscore pident qlen qstart qend slen sstart send" -max_target_seqs 1000 -evalue 10

