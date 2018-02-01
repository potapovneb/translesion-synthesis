#!/bin/bash

samples=$1
command=$2

known_commands=(
    consensus
    chimeric
    mutation
    mutation-bwa
    summary-A1
)

### define top project directory
PROJDIR=$PWD

cd $PROJDIR

while read line
do
    if [[ ! $line =~ "SampleID" ]] && [[ ! $line =~ "#" ]]
    then
	### sample details
	sampleId=`echo $line | cut -d, -f1`
	amplicon=`echo $line | cut -d, -f7`
	collectionPathUri=`echo $line | cut -d, -f9`

	### preview run info
	echo ""
	echo "sampleId=$sampleId"
	echo "reference=$reference"
	echo "collectionPathUri=$collectionPathUri"

	### output directory
	rundir=`printf "%s/samples/%05i" $PROJDIR $sampleId`
	rm -rf "$rundir"
	mkdir -p "$rundir"

	### run data analysis
    	qsub -v root="$PROJDIR",amplicon="plasmid-1",sampleId="$sampleId",collectionPathUri="$collectionPathUri",np="10" \
	     -N "ccs$sampleId" \
	     -o $rundir/$command.log \
	     -j yes \
	     $PROJDIR/scripts/ccs2-$command.sh
    fi
done < $samples
