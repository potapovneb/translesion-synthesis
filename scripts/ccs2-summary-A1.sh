#!/bin/bash

#$ -S /bin/bash
#$ -P longrun
#$ -pe smp 2

cat $JOB_SCRIPT >> $SGE_STDOUT_PATH

echo ""
echo "Started on $(date) on $(uname -n)"

module load p7zip

reference="$root/references/$amplicon/sequence/$amplicon.fasta"

### run directory
rundir=`printf "%s/samples/%05i/summary-bwa-pacbio/1FS" $root $sampleId`
mkdir -p "$rundir"
cd "$rundir"

### tally up mutations
echo ""
sampledir=`printf "%s/samples/%05i" $root $sampleId`
"$root"/bin/ccs2-summary-A1.pl \
    --mutation-dir "mutation-bwa-pacbio" \
    --np $np \
    --qv 93 \
    --mapq 60 \
    --lb 40 \
    --ub 40 \
    --rlen-min-a 50 \
    --rlen-max-a 50 \
    "$reference" "$sampledir" >summary-A1-$np.csv 2>summary-A1-$np.log
echo "Task 1 completed at $(date)"

echo ""
echo "Finished on $(date)"
