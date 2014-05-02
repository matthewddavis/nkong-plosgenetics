#!/usr/bin/python
'''
Written by Matt Davis - Feb 2012.
This script trims unzipped FASTQ files as a quick 
heurisitic for improving mapping quality with an alignment
program such as Bowtie or BWA.
'''
import sys


def trim_fastqfile(fastqfn, start, stop):
    '''Trims reads in a fastq file from start to stop. 
    This is used to get shorter (but still likely unique) reads 
    with high qscores (e.g. bases 10-30 for 20bp reads).

    Accepts:
    fastqfn - the name of the fastq file
    start - the base to start trimming at
    stop - the last base to include

    example - 
    trim_fastqfile('s1_sequence.fastq', 's1_trimmed.fastq', 10, 30)
 
    '''
    try:
        fastqfh = open(fastqfn, 'r')
    except:
        print "Error opening FASTQ file: ", fastqfn
        return

    try:
        outfn = fastqfn.split('.fastq')[0] + '.trimmed.fastq'
        outfh = open(outfn, 'w')
    except:
        print "Error opening output file: ", outfn
        sys.exit()

    print "Trimming reads from bp %i to %i for file %s..." % (start, stop, fastqfn)
    count = 0
    while (True):
        for i in range(4):
            seqname = fastqfh.readline().strip()
            if seqname == '': return
            seq = fastqfh.readline().strip()
            comment = fastqfh.readline().strip()
            qscore = fastqfh.readline().strip()
            outfh.write(seqname + '\n')
            outfh.write(seq[start:stop] + '\n')
            outfh.write(comment + '\n')
            outfh.write(qscore[start:stop] + '\n')
            count += 1
        if count % 1e6 == 0:
            print count, seq[start:stop]
    outfh.close()
    fastqfh.close()

def main():
    if len(sys.argv) < 2:
        print "Usage:  trimReads.py <infiles>"
        sys.exit()
    print sys.argv[1:]
    for fastqfn in sys.argv[1:]:
        start = 20
        stop = 50
        trim_fastqfile(fastqfn, start, stop)


if __name__ == '__main__':
    main()
