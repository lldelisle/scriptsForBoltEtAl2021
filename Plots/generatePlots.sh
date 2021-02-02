gitHubDirectory="/home/ldelisle/softwares/scriptsForBoltEtAl2021/"
pathWithCool="/scratch/ldelisle/BoltGEO/CHIC-seq/"
pathWithCR="/scratch/ldelisle/BoltGEO/CUTandRUN/"
pathWithATAC="/scratch/ldelisle/BoltGEO/ATAC-seq/"
plottingDirectory="/scratch/ldelisle/BoltPlot/"

echo "[annot]
file = ${gitHubDirectory}/scripts/annotations/HoxD_Elements_mm10_Minimal.bed
title = HoxD Elements
height = 1.5
fontsize = 10
file_type = bed
color = black
display = interleaved
labels = true
" > ${plottingDirectory}/annotations.ini

# Figure 2a
for tissue in DFL PFL; do
  ini_file=${plottingDirectory}/fig2a_${tissue}.ini
  echo "[hic]
file = ${pathWithCool}/E12.5_wt_${tissue}_CHIC.cool
title = ${tissue} 5kb
depth = 1300000
show_masked_bins = false
colormap = Greys
transform = log

[CnR]
file = ${pathWithCR}/E12.5_wt_${tissue}_CR_CTCF.bw
title = $tissue CTCF
height = 2
color = Black
min_value = 0
max_value = 120

[spacer]" > ${ini_file}
  cat ${plottingDirectory}/annotations.ini >> ${ini_file}
  pyGenomeTracks --tracks ${ini_file} \
    -o ${ini_file/ini/pdf} --region chr2:73805207-75678600 \
    --dpi 300 --fontSize 12
done

# Figure 2Sa

for tissue in DFL PFL FB; do
  ini_file=${plottingDirectory}/fig2Sa_${tissue}.ini
  echo "" > ${ini_file}
  for genotype in wt inv2; do
    echo "[hic ${genotype}]
file = ${pathWithCool}/E12.5_${genotype}_${tissue}_CHIC.cool
title = ${genotype} $tissue 5kb
depth = 750000
show_masked_bins = false
colormap = Greys
transform = log

[spacer]
" >> ${ini_file}
  done
  for genotype in wt inv2; do
    echo "[insulation ${genotype}]
file = ${plottingDirectory}/E12.5_${genotype}_${tissue}_CHIC.240kb_tad_score.bm
title = ${genotype} $tissue
type = line
file_type = bedgraph
height = 3
max_value = 1
" >> ${ini_file}
    if [ $genotype = "wt" ]; then
      echo "
color = blue
" >> ${ini_file}
    else
      echo "
color = red
overlay_previous = share-y
" >> ${ini_file}
    fi
  done
  echo "[spacer]
    
" >> ${ini_file}
  if [ $tissue = "FB" ]; then
    echo "[ATAC wt ${tissue}]
file = ${pathWithATAC}/E13.5_wt_${tissue}_ATAC.bw
title = ATAC wt ${tissue}
height = 3
color = Black
min_value=0
max_value=20

[spacer]

[ATAC inv2 DFL]
file = ${pathWithATAC}/E12.5_inv2_DFL_ATAC.bw
title = ATAC inv2 DFL
height = 3
color = Black
min_value = 0
max_value = 20
    " >> ${ini_file}
  else
    echo "[ATAC wt ${tissue}]
file = ${pathWithATAC}/E12.5_wt_${tissue}_ATAC.bw
title = ATAC wt ${tissue}
height = 3
color = Black
min_value=0
max_value=20

[spacer]

[ATAC inv2 ${tissue}]
file = ${pathWithATAC}/E12.5_inv2_${tissue}_ATAC.bw
title = ATAC inv2 ${tissue}
height = 3
color = Black
min_value = 0
max_value = 10
    " >> ${ini_file}
  fi
  echo "[spacer]

[CTCF wt DFL]
file = ${pathWithCR}/E12.5_wt_DFL_CR_CTCF.bw
title = CTCF wt DFL
height = 3
color = Black
min_value = 0
max_value = 120

[spacer]
" >> ${ini_file}
  cat ${plottingDirectory}/annotations.ini >> ${ini_file}
  pyGenomeTracks --tracks ${ini_file} \
    -o ${ini_file/ini/pdf} --region chr2:74277600-75147000 \
    --dpi 300 --fontSize 12
done

# Figure 2b
# Get public dataset:
wget "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2713nnn/GSM2713704/suppl/GSM2713704%5FPFL%5FE12%5FWt%5FH3K27ac%2EbedGraph%2Egz" -O ${plottingDirectory}/GSM2713704_PFL_E12_Wt_H3K27ac.bedGraph.gz
# Convert to bigwig
gunzip ${plottingDirectory}/GSM2713704_PFL_E12_Wt_H3K27ac.bedGraph.gz
cat ${plottingDirectory}/GSM2713704_PFL_E12_Wt_H3K27ac.bedGraph | awk '{print $0 > "tmp_"$1}'
cat tmp_* > ${plottingDirectory}/GSM2713704_PFL_E12_Wt_H3K27ac_sorted.bedGraph
rm tmp_*
bedGraphToBigWig ${plottingDirectory}/GSM2713704_PFL_E12_Wt_H3K27ac_sorted.bedGraph ~/genomes/fasta/mm10.fa.fai ${plottingDirectory}/GSM2713704_PFL_E12_Wt_H3K27ac.bw

ini_file=${plottingDirectory}/fig2b.ini
echo "[CnR]
file = ${pathWithCR}/E12.5_wt_PFL_CR_CTCF.bw
title = CTCF PFL
height = 1.5
color = Black
min_value = 0
max_value = 120


[spacer]

[Chip]
file = ${plottingDirectory}/GSM2713704_PFL_E12_Wt_H3K27ac.bw
title = H3K27ac PFL
height = 1.5
color = Black
min_value=0
max_value=100

[spacer]

[ATAC]
file = ${pathWithATAC}/E12.5_wt_PFL_ATAC.bw
title = ATAC wt PFL
height = 1.5
color = Black
min_value=0
max_value=20

[spacer]

[ATAC]
file = ${pathWithATAC}/E12.5_inv2_PFL_ATAC.bw
title = ATAC inv2 PFL
height = 1.5
color = Black
min_value=0
max_value=10

[spacer]
" > ${ini_file}
cat ${plottingDirectory}/annotations.ini >> ${ini_file}
pyGenomeTracks --tracks ${ini_file} \
  -o ${ini_file/ini/pdf} --region chr2:75129600-75677400 \
  --dpi 300 --fontSize 12

# Figure 2Sb
# Get public dataset:
wget "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM2713703&format=file&file=GSM2713703%5FDFL%5FE12%5FWt%5FH3K27ac%2EbedGraph%2Egz" -O ${plottingDirectory}/GSM2713704_PFL_E12_Wt_H3K27ac.bedGraph.gz
# Convert to bigwig
gunzip ${plottingDirectory}/GSM2713703_DFL_E12_Wt_H3K27ac.bedGraph.gz
cat ${plottingDirectory}/GSM2713703_DFL_E12_Wt_H3K27ac.bedGraph | awk '{print $0 > "tmp_"$1}'
cat tmp_* > ${plottingDirectory}/GSM2713703_DFL_E12_Wt_H3K27ac_sorted.bedGraph
rm tmp_*
bedGraphToBigWig ${plottingDirectory}/GSM2713703_DFL_E12_Wt_H3K27ac_sorted.bedGraph ~/genomes/fasta/mm10.fa.fai ${plottingDirectory}/GSM2713703_DFL_E12_Wt_H3K27ac.bw

ini_file=${plottingDirectory}/fig2Sb.ini
echo "[CnR]
file = ${pathWithCR}/E12.5_wt_DFL_CR_CTCF.bw
title = CTCF DFL
height = 1.5
color = Black
min_value = 0
max_value = 120


[spacer]

[Chip]
file = ${plottingDirectory}/GSM2713703_DFL_E12_Wt_H3K27ac.bw
title = H3K27ac DFL
height = 1.5
color = Black
min_value=0
max_value=100

[spacer]

[ATAC]
file = ${pathWithATAC}/E12.5_wt_DFL_ATAC.bw
title = ATAC wt DFL
height = 1.5
color = Black
min_value=0
max_value=20

[spacer]

[ATAC]
file = ${pathWithATAC}/E12.5_inv2_DFL_ATAC.bw
title = ATAC inv2 DFL
height = 1.5
color = Black
min_value=0
max_value=10

[spacer]
" > ${ini_file}
cat ${plottingDirectory}/annotations.ini >> ${ini_file}
pyGenomeTracks --tracks ${ini_file} \
  -o ${ini_file/ini/pdf} --region chr2:73900700-74703000 \
  --dpi 300 --fontSize 12
