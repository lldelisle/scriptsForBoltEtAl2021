# Using deepTools version 3.5 
pathWithBigWig="/scratch/ldelisle/BoltGEO/ATAC-seq/"
plottingDirectory="/scratch/ldelisle/BoltPlot/"

# Generate the matrix with the mean coverage in each ATAC peak
multiBigwigSummary BED-file \
  -b \
  "${pathWithBigWig}/E12.5_wt_PFL_ATAC.bw" \
  "${pathWithBigWig}/E12.5_inv2_PFL_ATAC.bw" \
  "${pathWithBigWig}/E12.5_wt_DFL_ATAC.bw" \
  "${pathWithBigWig}/E12.5_inv2_DFL_ATAC.bw" \
  "${pathWithBigWig}/E13.5_wt_FB_ATAC.bw" \
  -o "${plottingDirectory}/Figure2Sd_left.npz" \
  --outRawCounts 	"${plottingDirectory}/Figure2Sd_left_data.txt" \
  --numberOfProcessors 12 \
  --verbose \
  --BED "${plottingDirectory}/all_wt_ATAC_merged_HoxD_region_noGeneBodies.bed" \

# Generate the matrix with the mean coverage in the HoxD promoters
multiBigwigSummary BED-file \
  -b \
  "${pathWithBigWig}/E12.5_wt_PFL_ATAC.bw" \
  "${pathWithBigWig}/E12.5_inv2_PFL_ATAC.bw" \
  "${pathWithBigWig}/E12.5_wt_DFL_ATAC.bw" \
  "${pathWithBigWig}/E12.5_inv2_DFL_ATAC.bw" \
  "${pathWithBigWig}/E13.5_wt_FB_ATAC.bw" \
  -o "${plottingDirectory}/Figure2Sd_right.npz" \
  --outRawCounts 	"${plottingDirectory}/Figure2Sd_right_data.txt" \
  --numberOfProcessors 12 \
  --verbose \
  --BED "${plottingDirectory}/promoters_wt_ATAC_merged_Evx2-Hoxd1.bed" 

# R 3.6.0
Rscript ${gitHubDirectory}/ATACseq/05_figS2c_plot_clustering.R ${plottingDirectory}
