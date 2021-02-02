#!/bin/bash

#SBATCH -o slurm-%x-%A.out
#SBATCH -e slurm-%x-%A.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lucille.delisle@epfl.ch
#SBATCH --nodes 1
#SBATCH --mem 80G
#SBATCH --cpus-per-task 32
#SBATCH --time 48:00:00
#SBATCH --job-name cellranger3.1.0
#SBATCH --chdir /scratch/ldelisle/scRNA_3.1.0/

pathForCellRanger="~/cellranger/"

mkdir -p $pathForCellRanger
cd $pathForCellRanger

# Downlaod the software
if [ ! -e cellranger-3.1.0 ]; then
  wget -O cellranger-3.1.0.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-3.1.0.tar.gz?Expires=1611355081&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cHM6Ly9jZi4xMHhnZW5vbWljcy5jb20vcmVsZWFzZXMvY2VsbC1leHAvY2VsbHJhbmdlci0zLjEuMC50YXIuZ3oiLCJDb25kaXRpb24iOnsiRGF0ZUxlc3NUaGFuIjp7IkFXUzpFcG9jaFRpbWUiOjE2MTEzNTUwODF9fX1dfQ__&Signature=jWxKCIRHw~aSUYPUmDX9GkNhiddR6vFnt66a7vMZ~g3UzQ2WkGJ38XVf54yKjhBIKrOwldUByBLQsoz6s3Qv8Fmu6ssCqTaJFLpVVcGgD9KeEaKK4c~BXjvE6tEGELjErYwvgXv8ag~ejvvCLgxJLDA6o~grkApz~kJovWjP1ypi4lxNRgMECWblSgw-jmUZhMxm7LWyCNmmZGRabQpI5EG7uh0stGkrVGiZoHaq5zGkyt0evrjoqGeYtXebzsSt7i61FE2qUWLni3slqel~fJqQuyylfCYaKmFf2haDceBOa9IMDjpp~FRxE3p1-KvOvZpPJU5onsPye-19fybBuA__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA"
  tar -xzvf cellranger-3.1.0.tar.gz
fi

# Downlaod the gtf:
if [ ! -e mergeOverlapGenesOfFilteredTranscriptsOfMus_musculus.GRCm38.98_ExonsOnly_UCSC.gtf ]; then
  wget "https://zenodo.org/record/4456702/files/mergeOverlapGenesOfFilteredTranscriptsOfMus_musculus.GRCm38.98_ExonsOnly_UCSC.gtf.gz?download=1"
  gunzip mergeOverlapGenesOfFilteredTranscriptsOfMus_musculus.GRCm38.98_ExonsOnly_UCSC.gtf.gz
fi

module purge
export PATH=${pathForCellRanger}/cellranger-3.1.0/:$PATH

cellranger mkref --genome=merged.filtered.ensembl.grc38.98 \
  --fasta=/home/ldelisle/genomes/fasta/mm10.fa \
  --genes=mergeOverlapGenesOfFilteredTranscriptsOfMus_musculus.GRCm38.98_ExonsOnly_UCSC.gtf
