# Using deepTools version 3.5 
pathWithBigWig="/scratch/ldelisle/BoltGEO/CUTandRUN/"
plottingDirectory="/scratch/ldelisle/BoltPlot/"
gitHubDirectory="/home/softwares/scriptsForBoltEtAl2021/"

multiBigwigSummary bins \
  -b \
  "${pathWithBigWig}/E12.5_inv2_PPFL_CR_HOXD13_rep1.bw" \
  "${pathWithBigWig}/E12.5_inv2_PPFL_CR_HOXD13_rep2.bw" \
  "${pathWithBigWig}/E12.5_wt_DFL_CR_HOXD13_rep1.bw" \
  "${pathWithBigWig}/E12.5_wt_DFL_CR_HOXD13_rep2.bw" \
  "${pathWithBigWig}/E11.5_wt_WFL_ChIP_HOXA11_rep1.bw" \
  "${pathWithBigWig}/E11.5_wt_WFL_ChIP_HOXA11_rep2.bw" \
  -o "${plottingDirectory}/plotCorrelation_HOXD13_HOXA11_wt_v_inv2_allGenomeRegions.npz" \
  --numberOfProcessors 12 \
  --verbose

plotCorrelation --corData "${plottingDirectory}/plotCorrelation_HOXD13_HOXA11_wt_v_inv2_allGenomeRegions.npz" \
  --corMethod spearman \
  --skipZeros \
  --colorMap Blues \
  --plotNumbers \
  --labels "HOXD13 P-PFL Inv2" "HOXD13 P-PFL Inv2" "HOXD13 DFL wt" "HOXD13 DFL wt" "HOXA11 WFL wt" "HOXA11 WFL wt" \
  --whatToPlot "heatmap" \
  --plotHeight 6 \
  --plotWidth 7 \
  --plotFile "${plottingDirectory}/Figure4c.pdf"
