# This analysis was done on galaxy, here are the command lines used and the verison of the tools:
# HiCUP version 0.6.1
# Bowtie2 version 2.2.6
# Samtools version 1.2
hicup_digester --re1 '^GATC' --genome 'mm10' mm10.fa
mv *Digest_* digest_file.txt
hicup --zip --threads ${GALAXY_SLOTS:-1} --digest digest_file.txt --index 'bowtie2_index/mm10_UCSC/mm10_UCSC' --bowtie2 $BOWTIE_PATH_BASH --keep    dataset1.gz dataset2.gz 

# Custom script to get juicebox format (available in the scripts/hic_scripts/ directory):
python fromHicupToJuicebox.py --fragmentFile digest_file.txt --colForChr 1 --colForStart 2 --colForEnd 3 --colForID 4 --lineToSkipInFragmentFile 2 --methodForFrag hicup --useMid --output validPairs.txt output.bam
bash switchAndSort.sh validPairs.txt validPairs_sorted.txt

# Filter both valid pairs above mapq30:
# Version 1.1.0
python tools/stats/filtering.py validPairs_sorted.txt validPairs_sorted_filtered1.txt "c10__gt__=30 and c11__gt__=30" 11 "str,int,str,int,int,int,str,int,int,int,int" 0

# Filter both pairs in capture region
# Version 1.1.0
python tools/stats/filtering.py validPairs_sorted_filtered1.txt validPairs_sorted_filtered2.txt "(c3==__sq__chr2__sq__ and c4__lt__77000000 and c4__gt__72402000) and (c7==__dq__chr2__dq__ and c8__lt__77000000 and c8__gt__72402000)" 11 "str,int,str,int,int,int,str,int,int,int,int" 0

# cooler version 0.7.4
cooler csort -i tabix -c1 3 -c2 7 -p1 4 -p2 8 -o validPairs.csort.gz validPairs_sorted_filtered2.txt mm10.size

# cooler version 0.7.4
cooler makebins -o bins.5kb mm10.size 5000

cooler cload tabix --assembly mm10 -c2 7 -p2 8 bins.5kb validPairs.csort.gz output.cool

cooler balance --mad-max 5 --min-nnz 10 --min-count 0 --ignore-diags 2 --tol 1e-05 --max-iters 200 --cis-only -f output.cool

# Custom script to export to txt (available in the scripts/hic_scripts/ directory):
python cooler_getMatrix.py --input output.cool --output output.txt --balance --header --r1 chr2:72000000-77000000
