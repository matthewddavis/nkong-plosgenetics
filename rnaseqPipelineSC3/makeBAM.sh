#!/bin/bash

SAMPLEDIR=$1
READTYPE=$2

if [ $READTYPE == 'single' ]; then
    BWADIR='bwaSE'
elif [ $READTYPE == 'paired' ]; then
    BWADIR='bwaPE'
else
    echo -e "Usage:  <read type> should be either 'single or 'paired'"
    exit
fi

samfn=`echo $SAMPLEDIR/$BWADIR/*sam`
#note that the bamfn is the prefix without the .bam due to 
# how samtools sort works...
bamfn=`echo $samfn | sed 's/.sam//g'`
echo -e 'Converting' $samfn 'to .bam format...'
samtools view -bS $samfn > $bamfn.bam
echo -e '\t...done.'
    
echo -e 'Sorting' $bamfn.bam '...'
samtools sort $bamfn.bam $bamfn
echo -e '\t...done.'

