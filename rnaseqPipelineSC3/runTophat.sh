#!/bin/bash

SAMPLEDIR=$1

TOPHATDIR=$SAMPLEDIR/'tophat'
FASTQFILES=`ls $SAMPLEDIR/*.fastq`
echo $FASTQFILES

TOPHAT=`which tophat`
BOWTIEIDX='../reference/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome'
GTFFILE='../reference/Mus_musculus/UCSC/mm10/Annotation/Archives/archive-2013-03-06-15-06-02/Genes/genes.gtf'

$TOPHAT -p8 --no-novel-juncs --segment-length=18 -o $TOPHATDIR -G $GTFFILE $BOWTIEIDX $FASTQFILES
