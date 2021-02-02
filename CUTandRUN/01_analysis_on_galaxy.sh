# This analysis was done on galaxy, here are the command lines used and the verison of the tools:

# For the CUT&RUN:
# cutadapt version 1.16
cutadapt  -j ${GALAXY_SLOTS:-1}     -a 'Truseq R1'='GATCGGAAGAGCACACGTCTGAACTCCAGTCAC'         -A 'TruSeq R2'='GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT'      --output='out1.gz' --paired-output='out2.gz'  --error-rate=0.1 --times=1 --overlap=3         --minimum-length=15 --pair-filter=any   --quality-cutoff=30   'R1.fastq.gz' 'R2.fastq.gz'  > report.txt

# bowtie2  version 2.3.4.1
bowtie2  -p ${GALAXY_SLOTS:-4}  -x 'mm10_UCSC'   -1 'input_f.fastq.gz' -2 'input_r.fastq.gz' -I 0 -X 1000 --fr --no-mixed --no-discordant --dovetail                --very-sensitive  2> 'mapping_stat'  | samtools sort -@${GALAXY_SLOTS:-2} -O bam -o 'output.bam'

# samtools version 1.2
samtools view -o 'filtered.bam' -h   -b  -q 30 -f 0x2 input.bam 2>&1

# Picards version 1.56.0:
python devteam/picard/bf1c3f9f8282/picard/picard_wrapper.py -i "filtered.bam" -n "Dupes Marked " --tmpdir "/tmp" -o "rmdup.bam" --remdups "true" --assumesorted "true" --readregex "[a-zA-Z0-9]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*" --optdupdist "100" -j "$JAVA_JAR_PATH/MarkDuplicates.jar" -d "output" -t "output.txt" -e "bam"

# Bedtools version 2.18,2
bedtools bamtobed  -i ./input.bam > output.bed

# macs2 version 2.1.1.20160309
macs2 callpeak  --name 'MACS2'   -t 'output.bed'    --format BED   --gsize '1870000000'          --call-summits  --keep-dup '1'   --bdg  --qvalue '0.05'  --nomodel --extsize '200' --shift '-100'  2>&1 > macs2_stderr

# Convert to bigwig wig_to_bigWig tool version 1.1.0
grep -v "^track" macs2.bedgraph | wigToBigWig stdin mm10.len output.bw -clip 2>&1 || echo "Error running wigToBigWig." >&2

# For the ChIP:
# seqtk version 1.3.0
seqtk sample -s 4 'original.fastq.gz' 25000000.0  | gzip --no-name > 'fastq.gz'

# cutadapt version 1.16
cutadapt  -j ${GALAXY_SLOTS:-1}     -a 'Truseq R1'='GATCGGAAGAGCACACGTCTGAACTCCAGTCAC'         -A 'TruSeq R2'='GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT'      --output='out1.gz' --paired-output='out2.gz'  --error-rate=0.1 --times=1 --overlap=3         --minimum-length=15 --pair-filter=any   --quality-cutoff=30   'R1.fastq.gz' 'R2.fastq.gz'  > report.txt

# bowtie2 version 2.3.4.1
bowtie2  -p ${GALAXY_SLOTS:-4}  -x 'mm10_UCSC'   -1 'input_f.fastq.gz' -2 'input_r.fastq.gz'               2> 'mapping_stat'  | samtools sort -@${GALAXY_SLOTS:-2} -O bam -o 'output.bam'

# samtools version 1.2
samtools view -o 'filtered.bam' -h   -b  -q 30 -f 0x2 input.bam

# macs2 version 2.1.1.20160309
macs2 callpeak  --name 'MACS2'   -t 'filtered.bam'    --format BAMPE   --gsize '1870000000'          --call-summits  --keep-dup '1'   --bdg  --qvalue '0.05'  --mfold '5' '50'  --bw '300'  2>&1 > macs2_stderr

# Convert to bigwig wig_to_bigWig tool version 1.1.0
grep -v "^track" macs2.bedgraph | wigToBigWig stdin mm10.len output.bw -clip 2>&1 || echo "Error running wigToBigWig." >&2
