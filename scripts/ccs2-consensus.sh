#!/bin/bash

#$ -S /bin/bash
#$ -P longrun
#$ -pe smp 8

cat $JOB_SCRIPT >> $SGE_STDOUT_PATH

echo ""
echo "Started on $(date) on $(uname -n)"

module load smrtlink-3.1.1
module load samtools

reference="$root/references/$amplicon/sequence/$amplicon.fasta"

### create run directory
rundir=`printf "%s/samples/%05i" $root $sampleId`
echo $rundir
mkdir -p "$rundir"

### switch to working directory
cd "$rundir"

### look up sequencing data
echo ""
find "$collectionPathUri" -name "*.bax.h5" | sort -u > input.fofn
echo "Task 1 completed at $(date)"

### convert
echo ""
bax2bam --fofn=input.fofn -o movie
echo "Task 2 completed at $(date)"

### build ccs
echo ""
ccs --reportFile=subreads_ccs.csv --logFile=subreads_ccs.log --numThreads=7 --minPasses=1 movie.subreads.bam subreads_ccs.bam
echo "Task 3 ($strand) completed at $(date)"

### cleanup
echo ""
rm -f movie.*
rm -f subreads.*
echo "Task 4 completed at $(date)"

echo ""
echo "Finished on $(date)"
