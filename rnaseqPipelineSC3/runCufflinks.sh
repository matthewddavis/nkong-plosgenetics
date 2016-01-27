#!/bin/bash

SAMPLEDIR=$1

CUFFLINKS=`which cufflinks`
INPUTFILE=$SAMPLEDIR/tophat/accepted_hits.bam
OUTDIR=$SAMPLEDIR/cufflinks

#REFGTFFILE='/home/matt/data/Reference/IlluminaReferences/dm3/Drosophila_melanogaster/UCSC/dm3/Annotation/Genes/genes.gtf'
REFGTFFILE='../reference/Mus_musculus/UCSC/mm10/Annotation/Archives/archive-2013-03-06-15-06-02/Genes/genes.gtf'

# warning -- uses 7 cpu cores
$CUFFLINKS -p8 -G $REFGTFFILE -o $OUTDIR $INPUTFILE 



