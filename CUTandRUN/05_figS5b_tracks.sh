gitHubDirectory="/home/ldelisle/softwares/scriptsForBoltEtAl2021/"
pathWithCR="/scratch/ldelisle/BoltGEO/CUTandRUN/"
plottingDirectory="/scratch/ldelisle/BoltPlot/"

samples=('E12.5_wt_DFL_CR_HOXD13' 'E12.5_inv2_PPFL_CR_HOXD13' 'E11.5_wt_WFL_ChIP_HOXA11')
colors=('purple' 'grey' 'blue')
ini_file=${plottingDirectory}/figS5b.ini
echo "" > ${ini_file}
for index in ${!samples[@]}; do
  sample=${samples[${index}]}
  color=${colors[${index}]}
  for rep in 1 2; do
    echo "[${index}_rep${rep}_a]
file = ${pathWithCR}/${sample}_rep${rep}.bw
title = ${sample//_/ } rep${rep}
height = 3
color = ${color}
min_value=0
max_value=100

[spacer]

[${index}_rep${rep}_b]
file = ${pathWithCR}/${sample}_rep${rep}.narrowPeak
title = ${sample//_/ } rep${rep}
height = 1
color = ${color}
show_labels = false
type = box
file_type = narrow_peak

[spacer]
" >> ${ini_file}
  done
done
pyGenomeTracks --tracks ${ini_file} -o ${ini_file/.ini/.pdf} \
  --BED ${gitHubDirectory}/CUTandRUN/regions_to_plot.bed --dpi 300 \
  --fontSize 12 --plotWidth 15
