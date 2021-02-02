# Using deepTools version 3.5 
pathWithBigWig="/scratch/ldelisle/BoltGEO/ATAC-seq/"
plottingDirectory="/scratch/ldelisle/BoltPlot/"

# Figure 2c left: ATAC peaks around HoxD but not in gene bodies:
computeMatrix reference-point \
  --referencePoint center \
  -b 1000 -a 1000 --verbose -p 12 \
  --samplesLabel 'PFL wt' 'PFL inv2' 'DFL wt' 'DFL inv2' 'FB wt' \
  --sortRegions keep \
  -S \
  "${pathWithBigWig}/E12.5_wt_PFL_ATAC.bw" \
  "${pathWithBigWig}/E12.5_inv2_PFL_ATAC.bw" \
  "${pathWithBigWig}/E12.5_wt_DFL_ATAC.bw" \
  "${pathWithBigWig}/E12.5_inv2_DFL_ATAC.bw" \
  "${pathWithBigWig}/E13.5_wt_FB_ATAC.bw" \
  -R \
  "${plottingDirectory}/all_wt_ATAC_merged_HoxD_region_noGeneBodies.bed" \
  -o "${plottingDirectory}/all_wt_ATAC_merged_HoxD_region_noGeneBodies_centered1000_matrix"

plotHeatmap --colorMap Blues Blues Reds Reds Greys \
  --missingDataColor 1 \
  --sortRegions no \
  --zMin 0 0 0 0 0 --zMax 40 20 60 20 40 \
  -m "${plottingDirectory}/all_wt_ATAC_merged_HoxD_region_noGeneBodies_centered1000_matrix" \
  --plotType lines \
  --heatmapHeight 35 \
  --outFileName "${plottingDirectory}/Figure2c_left.pdf"

# Figure 2c right: promoters in HoxD cluster with ATAC peak
computeMatrix reference-point \
  --referencePoint center \
  -b 1000 -a 1000 --verbose -p 12 \
  --samplesLabel 'PFL wt' 'PFL inv2' 'DFL wt' 'DFL inv2' 'FB wt' \
  --sortRegions keep \
  -S \
  "${pathWithBigWig}/E12.5_wt_PFL_ATAC.bw" \
  "${pathWithBigWig}/E12.5_inv2_PFL_ATAC.bw" \
  "${pathWithBigWig}/E12.5_wt_DFL_ATAC.bw" \
  "${pathWithBigWig}/E12.5_inv2_DFL_ATAC.bw" \
  "${pathWithBigWig}/E13.5_wt_FB_ATAC.bw" \
  -R \
  "${plottingDirectory}/promoters_wt_ATAC_merged_Evx2-Hoxd1.bed" \
  -o "${plottingDirectory}/promoters_wt_ATAC_merged_Evx2-Hoxd1_centered1000_matrix"

plotHeatmap --colorMap Blues Blues Reds Reds Greys \
  --missingDataColor 1 \
  --sortRegions no \
  --zMin 0 0 0 0 0 --zMax 40 20 60 20 40 \
  -m "${plottingDirectory}/promoters_wt_ATAC_merged_Evx2-Hoxd1_centered1000_matrix" \
  --plotType lines \
  --heatmapHeight 35 \
  --outFileName "${plottingDirectory}/Figure2c_right.pdf"
