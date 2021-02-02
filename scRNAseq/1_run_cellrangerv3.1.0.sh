#!/bin/bash

#SBATCH -o slurm-%x-%A_%2a.out
#SBATCH -e slurm-%x-%A_%2a.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lucille.delisle@epfl.ch
#SBATCH --nodes 1
#SBATCH --mem 80G
#SBATCH --cpus-per-task 32
#SBATCH --time 48:00:00
#SBATCH --array=1-3
#SBATCH --job-name cellranger3.1.0
#SBATCH --chdir /scratch/ldelisle/scRNA_3.1.0/

gitHubDirectory=$1
pathForCellRanger="~/cellranger/"

pathForTableWithSamples="${gitHubDirectory}/scripts/table_scRNAseq.txt" #PUT YOUR TABLE with one sample per line 
pathForFastq=$PWD/fastq/

sample=`cat $pathForTableWithSamples | awk -v i=$SLURM_ARRAY_TASK_ID 'NR==i{print $1}'`

module purge
export PATH=${pathForCellRanger}/cellranger-3.1.0/:$PATH

cellranger count --id=$sample \
                   --transcriptome=${pathForCellRanger}/merged.filtered.ensembl.grc38.98 \
                   --localcores=32 \
                   --fastqs=$pathForFastq/$sample \
                   --expect-cells=10000
