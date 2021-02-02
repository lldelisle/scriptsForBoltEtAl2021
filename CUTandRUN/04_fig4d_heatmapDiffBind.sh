gitHubDirectory="/home/ldelisle/softwares/scriptsForBoltEtAl2021/"
pathWithBigWig="/scratch/ldelisle/BoltGEO/CUTandRUN/"
plottingDirectory="/scratch/ldelisle/BoltPlot/"
# Using R version 3.6
Rscript ${gitHubDirectory}/scripts/runDiffBind.R ${gitHubDirectory}/scripts/tables/CUTandRUN_peaksets_A11PFL_D13wtDFL.csv ${plottingDirectory}/A11vsD13

# Using deepTools version 3.5
computeMatrix reference-point \
  --referencePoint center \
  -b 1000 -a 1000 --verbose -p 12 \
  --samplesLabel 'HOXA11 PFL wt' 'HOXD13 PFL inv2' 'HOXD13 DFL wt' \
  --sortRegions descend \
  --outFileSortedRegions "${plottingDirectory}/Figure2d_regionsOutput.bed" \
  -S \
  "${pathWithBigWig}/E11.5_wt_WFL_ChIP_HOXA11_rep2.bw" \
  "${pathWithBigWig}/E12.5_inv2_PPFL_CR_HOXD13_rep2.bw" \
  "${pathWithBigWig}/E12.5_wt_DFL_CR_HOXD13_rep2.bw" \
  -R \
  "${plottingDirectory}/A11vsD13_WFL.bed" \
  "${plottingDirectory}/A11vsD13_NDB.bed" \
  "${plottingDirectory}/A11vsD13_DFL.bed" \
  -o "${plottingDirectory}/CUTandRUN_peaksets_A11PFL_D13wtDFL_centered1000_matrix"

plotHeatmap --colorMap Blues Greys Purples \
  --missingDataColor 1 \
  -m "${plottingDirectory}/CUTandRUN_peaksets_A11PFL_D13wtDFL_centered1000_matrix" \
  --plotType lines \
  --outFileName "${plottingDirectory}/Figure4d.pdf"
