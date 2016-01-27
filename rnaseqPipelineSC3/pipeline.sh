#!/bin/bash

if [ $# -lt 1 ]; then
    echo -e 'Usage: ' `basename $0` '<sample directory>'
    exit
fi

SAMPLEDIR=$1

echo -e 'Mapping reads with Tophat for' $SAMPLEDIR'...\n'
./runTophat.sh $SAMPLEDIR

echo -e 'Making BigWig file for' $SAMPLEDIR'...\n'
./makeNormalizedBigWig.sh $SAMPLEDIR

echo -e "Running TopHat for" $SAMPLEDIR"...\n"
./runTophat.sh

echo -e "Running Cufflinks for" $SAMPLEDIR"...\n"
./runCufflinks.sh