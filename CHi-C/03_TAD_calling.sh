plottingDirectory="/scratch/ldelisle/BoltPlot/"
gitHubDirectory="/home/softwares/scriptsForBoltEtAl2021/"
pathWithCool="/scratch/ldelisle/BoltGEO/CHIC-seq/"

# HiCExplorer version 3.5.1
for f in ${pathWithCool}/*.cool; do
  hicFindTADs -m "$f" \
    --minBoundaryDistance 200000 \
    --correctForMultipleTesting fdr \
    --outPrefix ${plottingDirectory}/$(basename ${f/.cool/}).240kb \
    --minDepth 240000 \
    --maxDepth $((2 * 240000)) \
    --step $((2 * 240000)) \
    --numberOfProcessors 12
done
