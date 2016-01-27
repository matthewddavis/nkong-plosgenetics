#!/bin/bash
SAMPLEDIR=$1
CHROMSIZEFILE='../../reference/sacCer3.chrom.sizes'

sample=`basename $SAMPLEDIR`
OUTPUTDIR=$SAMPLEDIR/tophat
BAMFILE=$OUTPUTDIR/accepted_hits
EXTEND=200

echo -e 'Making .bai file for' $BAMFILE.bam '...'
samtools index $BAMFILE.bam $BAMFILE.bai
echo -e '\t...done.'

# Number of reads
echo 'Calculating librarysize...'
librarySize=$(samtools idxstats $OUTPUTDIR/accepted_hits.bam | awk '{total+=$3}END{print total}')
echo $librarySize > $OUTPUTDIR/librarySize
echo -e '\t done.'
echo -e 'Creating density.gz file...'
# Create density file: extend reads, calculate read density at each position and normalize the library size to 1 million reads
bamToBed -i $BAMFILE.bam | awk -vOFS='\t' '{print $1,$2,$3,$6}' | sort -k1,1 -k2,2n | genomeCoverageBed -i stdin -g $CHROMSIZEFILE -d | awk -vOFS='\t' -vSIZE=$librarySize '{print $1,$2,$2+1,$3*1000000/SIZE}' | gzip > $OUTPUTDIR/${sample}.density.gz
echo -e '\t done...'
# Create WIG file
echo -e 'Creating wig.gz file...'
gunzip -c $OUTPUTDIR/${sample}.density.gz | awk -vOFS='\t' '($4!=0) {if(!chrom[$1]){print "variableStep chrom="$1;chrom[$1]=1};print $2,$4}' | gzip > $OUTPUTDIR/${sample}.wig.gz
echo -e '\t done...'
# Create BigWig file
echo -e 'Creating BigWig file...'
./wigToBigWig $OUTPUTDIR/${sample}.wig.gz $CHROMSIZEFILE $OUTPUTDIR/${sample}.bw
echo -e '\tdone...'
echo -e 'Creating BedGraph file...'
./bigWigToBedGraph $OUTPUTDIR/${sample}.bw $OUTPUTDIR/${sample}.bedGraph
echo -e '\tdone...'

