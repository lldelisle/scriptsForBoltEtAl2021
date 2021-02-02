# Homer version 4.10
pathWithPeaks="/scratch/ldelisle/BoltGEO/CUTandRUN/"
plottingDirectory="/scratch/ldelisle/BoltPlot/"

for f in ${pathWithPeaks}/*HOXD13*_rep2.narrowPeak; do
  findMotifsGenome.pl \
    ${f} \
    mm10 \
    ${plottingDirectory}/$(basename ${f/.narrowPeak/_50bp}) \
    -size 50 -p 12 > ${plottingDirectory}/$(basename ${f/.narrowPeak/_50bp}).log 2>&1 &
done
