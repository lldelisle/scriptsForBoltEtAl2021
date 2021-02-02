# This analysis was done on galaxy, here are the command lines used and the verison of the tools:
# bowtie2  version 2.3.4.1
bowtie2  -p ${GALAXY_SLOTS:-4}  -x 'mm10_UCSC'   -1 'input_f.fastq.gz' -2 'input_r.fastq.gz' -I 0 -X 2000 --fr   --dovetail                --very-sensitive-local  2> 'mapping_stat'  | samtools sort -@${GALAXY_SLOTS:-2} -O bam -o 'file.bam'
# bamtools version 2.4.1
bamtools filter -script 'rules.txt' -in localbam.bam -out 'filtered.bam'
# Here is the content of rules.txt:
# {
#     "filters": [
#         {
#             "isProperPair": "true", 
#             "mapQuality": ">=30", 
#             "id": "1", 
#             "reference": "!chrM"
#         }
#     ]
# }
# Picards version 1.56.0:
python devteam/picard/bf1c3f9f8282/picard/picard_wrapper.py -i "filtered.bam" -n "Dupes Marked " --tmpdir "/tmp" -o "rmdup.bam" --remdups "true" --assumesorted "true" --readregex "[a-zA-Z0-9]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*" --optdupdist "100" -j "$JAVA_JAR_PATH/MarkDuplicates.jar" -d "output" -t "output.txt" -e "bam"
# Bedtools version 2.18.2
bedtools bamtobed  -i ./input.bam > output.bed

# Custom tool to get coverage (available in the scripts/hic_scripts/ directory):
python getTn5ExtendedCoverage.py --input input.bam --length 20 --output output.bedgraph

# Convert to bigwig wig_to_bigWig tool version 1.1.0
grep -v "^track" output.bedgraph | wigToBigWig stdin mm10.len output.bw -clip 2>&1 || echo "Error running wigToBigWig." >&2

# macs2 version 2.1.1.20160309
macs2 callpeak  --name 'MACS2'   -t 'output.bed'    --format BED   --gsize '2652783500'          --call-summits  --keep-dup 'all'    --qvalue '0.05'  --nomodel --extsize '200' --shift '-100'  2>&1 > macs2_stderr
