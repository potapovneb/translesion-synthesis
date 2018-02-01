#!/bin/bash

#$ -S /bin/bash
#$ -P longrun

cat $JOB_SCRIPT >> $SGE_STDOUT_PATH

echo ""
echo "Started on $(date) on $(uname -n)"

module load p7zip
module load bwa-0.7.13
module load samtools

reference="$root/references/$amplicon/sequence/$amplicon.fasta"

run directory
ccsdir=`printf "%s/samples/%05i" $root $sampleId`
rundir=`printf "%s/samples/%05i/mutation-bwa-pacbio" $root $sampleId`

### run directory
echo ""
mkdir -p "$rundir"
cd "$rundir"
echo "Task 1 completed at $(date)"

### copy and index reference fasta
echo ""
cp "$reference" reference.fasta
echo "Task 2 completed at $(date)"

### align reads
echo ""
"$root"/bin/bam2fq.pl "$ccsdir"/subreads_ccs.bam subreads_ccs.fastq
bwa index reference.fasta
bwa mem -x pacbio reference.fasta subreads_ccs.fastq > aligned_reads.sam
samtools view -Sb aligned_reads.sam -T reference.fasta > aligned_reads.bam
echo "Task 3 completed at $(date)"

### call mutations
echo ""
"$root"/bin/ccs2-mutations.pl aligned_reads.bam reference.fasta >variants.csv
echo "Task 4 completed at $(date)"

### pacbio read stats
echo ""
"$root"/bin/bam2csv.pl "$ccsdir"/subreads_ccs.bam zmws.csv
echo "Task 5 completed at $(date)"

### alignment stats
echo ""
"$root"/bin/alignment.pl aligned_reads.bam reference.fasta > aln.csv
echo "Task 6 completed at $(date)"

### compress
echo ""
7za a -t7z -m0=lzma -mx=9 -mfb=64 -md=32m -ms=on variants.csv.7z variants.csv
7za a -t7z -m0=lzma -mx=9 -mfb=64 -md=32m -ms=on zmws.csv.7z zmws.csv
7za a -t7z -m0=lzma -mx=9 -mfb=64 -md=32m -ms=on aln.csv.7z aln.csv
echo "Task 7 completed at $(date)"

### cleanup
echo ""
rm -f reference.*
rm -f variants.csv
rm -f zmws.csv
rm -f aln.csv
rm -f subreads_ccs.fastq
rm -f aligned_reads.sam
echo "Task 8 completed at $(date)"

echo ""
echo "Finished on $(date)"
