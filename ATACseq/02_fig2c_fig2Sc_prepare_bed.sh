# Using bedools version 2.27.1
pathWithPeaks="/scratch/ldelisle/BoltGEO/ATAC-seq/"
plottingDirectory="/scratch/ldelisle/BoltPlot/"
gitHubDirectory="/home/ldelisle/softwares/scriptsForBoltEtAl2021/"

# Merge all ATAC peaks from wt samples:
cat ${pathWithPeaks}/*wt*ATAC.narrowPeak | cut -f 1-3 | bedtools sort -i - > ${plottingDirectory}/all_wt_ATAC.bed
bedtools merge -i ${plottingDirectory}/all_wt_ATAC.bed > ${plottingDirectory}/all_wt_ATAC_merged.bed
# Remove random chromosomes or contigs
cat ${plottingDirectory}/all_wt_ATAC_merged.bed | grep -v random | grep -v chrUn > ${plottingDirectory}/all_wt_ATAC_merged_main.bed 
# 177568 regions

# Overlap with the HoxD region, approximately the 2 TADs:
echo -e "chr2\t73950000\t75655000" > ${plottingDirectory}/HoxD_region.bed
bedtools intersect -a ${plottingDirectory}/all_wt_ATAC_merged.bed \
  -b ${plottingDirectory}/HoxD_region.bed > ${plottingDirectory}/all_wt_ATAC_merged_HoxD_region.bed

# Downlaod the gtf:
wget https://zenodo.org/record/4456702/files/mergeOverlapGenesOfFilteredTranscriptsOfMus_musculus.GRCm38.98_ExonsOnly_UCSC.gtf.gz?download=1 -O ${plottingDirectory}/mergeOverlapGenesOfFilteredTranscriptsOfMus_musculus.GRCm38.98_ExonsOnly_UCSC.gtf.gz

# Reformat it:
zcat ${plottingDirectory}/mergeOverlapGenesOfFilteredTranscriptsOfMus_musculus.GRCm38.98_ExonsOnly_UCSC.gtf.gz \
  | sed -e 's/\"\"/\"/g' -e 's/\t\"/\t/g' -e 's/;\"/;/g' > ${plottingDirectory}/mergeOverlapGenesOfFilteredTranscriptsOfMus_musculus.GRCm38.98_ExonsOnly_UCSC_formatted.gtf

# Generate a bed file with gene bodies to exclude:
echo "track name=genes_to_exclude" > ${gitHubDirectory}/scripts/annotations/geneBodies.bed
cat ${plottingDirectory}/mergeOverlapGenesOfFilteredTranscriptsOfMus_musculus.GRCm38.98_ExonsOnly_UCSC_formatted.gtf | awk 'BEGIN{start=10000000000;end=0}$0~/\"Lnpk\"/{if($4<start){start=$4};if($5>end){end=$5};chrom=$1}END{print chrom"\t"start - 1"\t"end"\tLnpk"}' >> ${gitHubDirectory}/scripts/annotations/geneBodies.bed
cat ${plottingDirectory}/mergeOverlapGenesOfFilteredTranscriptsOfMus_musculus.GRCm38.98_ExonsOnly_UCSC_formatted.gtf | awk 'BEGIN{start=10000000000;end=0}$0~/\"Hoxd1\"/ || $0~/\"Evx2\"/{if($4<start){start=$4};if($5>end){end=$5};chrom=$1}END{print chrom"\t"start - 1"\t"end"\tEvx2-Hoxd1"}' >> ${gitHubDirectory}/scripts/annotations/geneBodies.bed
cat ${plottingDirectory}/mergeOverlapGenesOfFilteredTranscriptsOfMus_musculus.GRCm38.98_ExonsOnly_UCSC_formatted.gtf | awk 'BEGIN{start=10000000000;end=0}$0~/\"Mtx2\"/{if($4<start){start=$4};if($5>end){end=$5};chrom=$1}END{print chrom"\t"start - 1"\t"end"\tMtx2"}' >> ${gitHubDirectory}/scripts/annotations/geneBodies.bed


# Exclude the gene bodies
bedtools intersect -a ${plottingDirectory}/all_wt_ATAC_merged_HoxD_region.bed \
  -b ${gitHubDirectory}/scripts/annotations/geneBodies.bed -v > ${plottingDirectory}/all_wt_ATAC_merged_HoxD_region_noGeneBodies.bed

# Also generate a bed file with promoters (-1kb, +100bp from TSS)
cat ${plottingDirectory}/mergeOverlapGenesOfFilteredTranscriptsOfMus_musculus.GRCm38.98_ExonsOnly_UCSC_formatted.gtf | awk '$3 == "exon" && $0~/exon_number \"1\"/{if($7 == "+"){print $1"\t"$4-1"\t"$4"\t.\t0\t"$7}else{print $1"\t"$5-1"\t"$5"\t.\t0\t"$7}}' > ${plottingDirectory}/tss_with_strand.bed
bedtools slop -i ${plottingDirectory}/tss_with_strand.bed -g ~/genomes/fasta/mm10.fa.fai -l 1000 -r 100 -s > ${plottingDirectory}/promoters.bed

# Select the promoters which are part of ATAC_wt peaks and between Evx2 and Hoxd1:
bedtools intersect -a ${plottingDirectory}/promoters.bed -b ${plottingDirectory}/all_wt_ATAC.bed -wa -u > ${plottingDirectory}/promoters_wt_ATAC.bed
bedtools sort -i ${plottingDirectory}/promoters_wt_ATAC.bed | bedtools merge -i - > ${plottingDirectory}/promoters_wt_ATAC_merged.bed
grep "Evx2-Hoxd1" ${gitHubDirectory}/scripts/annotations/geneBodies.bed > ${plottingDirectory}/Evx2-Hoxd1.bed
bedtools intersect -a ${plottingDirectory}/promoters_wt_ATAC_merged.bed -b ${plottingDirectory}/Evx2-Hoxd1.bed > ${plottingDirectory}/promoters_wt_ATAC_merged_Evx2-Hoxd1.bed
