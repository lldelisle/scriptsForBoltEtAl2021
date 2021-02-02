# Using R version 3.6.3
plottingDirectory="/scratch/ldelisle/BoltPlot/"
gitHubDirectory="/home/ldelisle/softwares/scriptsForBoltEtAl2021/"
pathWithCool="/scratch/ldelisle/BoltGEO/CHIC-seq/"

# From the cool file:
for f in ${pathWithCool}/*.cool; do
  python ${gitHubDirectory}/scripts/hic_scripts/cooler_getMatrix.py --input $f --output $(basename ${f/cool/txt}) --balance --header --r1 chr2:72000000-77000000
done

# Plot the differential heatmaps
Rscript ${gitHubDirectory}/scripts/hic_scripts/plotHiC.R ${gitHubDirectory}/scripts/hic_scripts/Figure_2b_configFile_wtPL_against_invPL_wtMapping_Hoxd_against_CS39_Hnrnpa3_inclusive.R
Rscript ${gitHubDirectory}/scripts/hic_scripts/plotHiC.R ${gitHubDirectory}/scripts/hic_scripts/Figure_2Sb_configFile_wtDL_InvDL_subtraction_wtCDOM_region.R
