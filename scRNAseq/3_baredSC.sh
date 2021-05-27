gitHubDirectory="/home/ldelisle/softwares/scriptsForBoltEtAl2021/"
mkdir -p baredSC_2d
gunzip -c ${gitHubDirectory}/scRNAseq/meta_data_scRNAseq.txt.gz > meta_data_scRNAseq.txt
# Split by cluster
cat meta_data_scRNAseq.txt | awk -v output="cluster" -v OFS="\t" '
NR==1{
  for(k=1;k<=NF;k++){
    if($k == "seurat_clusters"){
      # The table was exported with row names:
      i=k+1
    }
  }
  header=$0
}
NR>1{
  cluster=$i
  if(! (cluster in seen)){
    # Print the header
    print header > output"_"cluster".txt"
    seen[cluster] = 1
  }
  # Remove the row name
  $1=""
  $0=$0
  $1=$1
  print $0 > output"_"cluster".txt"
}'
# Then we generate the table for parallel mcmc in 2d:
pathForTable="${gitHubDirectory}/scRNAseq/table_baredSC_2d.txt"
echo -e "cluster_3.txt\tHoxa11\tHoxd13\t3\t3\tgenotype" > $pathForTable
echo -e "cluster_3.txt\tHoxa13\tHoxd13\t3\t3\tgenotype" >> $pathForTable
echo -e "cluster_6.txt\tHoxa11\tHoxd13\t3\t3\tgenotype" >> $pathForTable
echo -e "cluster_6.txt\tHoxd11\tHoxd13\t3\t3\tgenotype" >> $pathForTable
# I launch the parallel mcmc:
sbatch --array 1-4 --chdir $PWD/ ${gitHubDirectory}/scRNAseq/sbatch_baredSC_2d.sh ${pathForTable} $PWD/baredSC_2d/mcmc/
# When it is finished I plot
Rscript ${gitHubDirectory}/scRNAseq/MCMC_2d_plots.R $PWD/ baredSC_2d/mcmc/ ${pathForTable} Figure3f

# Then we generate the table for parallel mcmc in 1d:
pathForTable="${gitHubDirectory}/scRNAseq/table_baredSC_1d.txt"
if [ -e $pathForTable ]; then
  rm $pathForTable
fi
nsim=0
for cluster in 1 6; do
  for gene in Hoxd13 Hoxa11 Hoxd11 Hoxa13; do
    echo -e "cluster_${cluster}.txt\t${gene}\t3\tgenotype" >> $pathForTable
    nsim=$((nsim + 1))
  done
done
# I launch the parallel mcmc:
sbatch --array 1-${nsim} --chdir $PWD/ ${gitHubDirectory}/scRNAseq/sbatch_baredSC_1d.sh ${pathForTable} $PWD/baredSC_1d/mcmc/
# When it is finished I plot
Rscript ${gitHubDirectory}/scRNAseq/MCMC_1d_plots.R $PWD/ baredSC_1d/mcmc/ ${pathForTable} Figure4Sb
