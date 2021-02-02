pathWithPeaks="/scratch/ldelisle/BoltGEO/CUTandRUN/"
plottingDirectory="/scratch/ldelisle/BoltPlot/"
gitHubDirectory="/home/ldelisle/softwares/scriptsForBoltEtAl2021/"

# bedtools 2.27.1
all_peaks=""
for cond in E11.5_wt_WFL_ChIP_HOXA11 E12.5_inv2_PPFL_CR_HOXD13 E12.5_wt_DFL_CR_HOXD13; do
  bedtools intersect -a ${pathWithPeaks}/${cond}_rep1.narrowPeak \
    -b ${pathWithPeaks}/${cond}_rep2.narrowPeak \
    | bedtools sort | cut -f 1-3 | uniq > ${plottingDirectory}/${cond}.bed
  all_peaks="$all_peaks ${plottingDirectory}/${cond}.bed"
done

wc -l $all_peaks
#   16740 /scratch/ldelisle/BoltPlot//E11.5_wt_WFL_ChIP_HOXA11.bed
#    2099 /scratch/ldelisle/BoltPlot//E12.5_inv2_PPFL_CR_HOXD13.bed
#   24141 /scratch/ldelisle/BoltPlot//E12.5_wt_DFL_CR_HOXD13.bed

for peak1 in $all_peaks; do
  for peak2 in $all_peaks; do
    if [ $peak1 = $peak2 ]; then
      continue
    fi
    output=${peak1}vs$(basename $peak2)
    if [ -e ${peak2}vs$(basename $peak1) ]; then
      continue
    fi
    bedtools intersect -a ${peak1} \
    -b ${peak2} > $output
  done
done
wc -l ${plottingDirectory}/*vs*.bed
#   1453 /scratch/ldelisle/BoltPlot//E11.5_wt_WFL_ChIP_HOXA11.bedvsE12.5_inv2_PPFL_CR_HOXD13.bed
#   6182 /scratch/ldelisle/BoltPlot//E11.5_wt_WFL_ChIP_HOXA11.bedvsE12.5_wt_DFL_CR_HOXD13.bed
#   1955 /scratch/ldelisle/BoltPlot//E12.5_inv2_PPFL_CR_HOXD13.bedvsE12.5_wt_DFL_CR_HOXD13.bed

bedtools multiinter -i ${all_peaks} | awk -v OFS="\t" '$4==3{print $1,$2,$3}' > ${plottingDirectory}/intersect3HOX.bed
wc -l ${plottingDirectory}/intersect3HOX.bed
# 1359 /scratch/ldelisle/BoltPlot//intersect3HOX.bed

# Export the number of peaks in txt file:
if [ ! -e ${plottingDirectory}/intersect_counts.txt ]; then
  wc -l $all_peaks > ${plottingDirectory}/intersect_counts.txt
  wc -l ${plottingDirectory}/*vs*.bed >> ${plottingDirectory}/intersect_counts.txt
  wc -l ${plottingDirectory}/intersect3HOX.bed >> ${plottingDirectory}/intersect_counts.txt
fi

Rscript ${gitHubDirectory}/scripts/plotEuler.R ${plottingDirectory}/intersect_counts.txt ${plottingDirectory}/Figure4e.pdf

