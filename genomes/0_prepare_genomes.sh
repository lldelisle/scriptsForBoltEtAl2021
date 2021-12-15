#!/bin/bash

#SBATCH -o slurm-%x-%A.out
#SBATCH -e slurm-%x-%A.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lucille.delisle@epfl.ch
#SBATCH --nodes 1
#SBATCH --mem 80G
#SBATCH --cpus-per-task 1
#SBATCH --time 48:00:00
#SBATCH --job-name preparegenome
#SBATCH --chdir /home/ldelisle/genomes/fasta/

# Get mm10 from UCSC:
wget http://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz
gunzip mm10.fa.gz

# To prepare the mutant genome creation
# Get the mutant chr2
wget "https://zenodo.org/record/4456654/files/chr2_InvCS65-SB7.fa.gz?download=1" -O chr2_inv2.fa.gz
gunzip chr2_inv2.fa.gz
# Get all chrs from mm10 except chr2 using seqtk
cat mm10.fa | grep ">" | grep -v chr2 | sed 's/>//' > listOfChrs.txt
# I reformat them:
# This can be long:
seqtk seq -U mm10.fa | seqtk subseq -l 60 - listOfChrs.txt > allChrsIncludingContigsExceptchr2.fa
cat chr2_inv2.fa allChrsIncludingContigsExceptchr2.fa > mm10_inv2.fa
