#!/usr/bin/python
import sys
from numpy import array, mean, zeros

def getChromSizes(chromfn):
    '''
    '''
    try:
        print "Opening chrom sizes file..."
        chromfh = open(chromfn)
    except:
        print "Can't open ", chromfn
    sizes = {}
    for line in chromfh:
        line = line.split()
        print line
        if line == []:
            continue
        sizes[line[0]] = int(line[1])
    return sizes

def adjustChromSizes(chromSizes, binsize):
    '''
    '''
    return

def binDensity(coveragefn, binsize, chromSizes):
    '''
    '''
    coveragefh = open(coveragefn)
    
    chromBins = {}
    for chrom in chromSizes:
        #chromBins[chrom] = array(zeros(chromSizes[chrom] / binsize))
        chromBins[chrom] = []
    #print chromBins
    while True:
        scores = []
        for i in range(binsize):
            line = coveragefh.readline().split()
            if line == []:
                scoremean = mean(scores)
                chromBins[chrom].append(scoremean)
                break
            chrom = line[0]
            if chrom not in chromBins:
                print "Error:  Chromosome in density file that is not in chromsizes file. exiting..."
                sys.exit()
            scores.append(float(line[-1]))
        scoremean = mean(scores)
        chromBins[chrom].append(scoremean)
        if line == []:
            break
    return chromBins

def writeWig(chromBins, chromSizes, binsize, wigfn):
    '''
    '''
    wigfh = open(wigfn, 'w')
    wigfh.write('track name=%s description=%s type=wiggle_0 visibility=full autoScale=off viewLimits=0.0:20.0 color=100,150,0 maxHeightPixels=50:50:11 yLineMark=$val yLineOnOff=on smoothingWindow=10\n' % ("testwig", "testDescription",))
    for chrom in sorted(chromBins):
        start = 1
        wigfh.write('fixedStep chrom=%s start=%s step=%s span=%s\n' % (chrom, start, binsize, binsize))
        for binval in chromBins[chrom]:
            if (start + binsize) > chromSizes[chrom]:
                print chrom, "is too long!", start, binsize
                continue
            wigfh.write('%s\n' % (binval))
            start += binsize

def main():
    coveragefn = sys.argv[1]
    wigfn = sys.argv[2]
    chromfn = 'mm10.chrom.sizes'
    binsize = 100
    chromSizes = getChromSizes(chromfn)
    chromBins = binDensity(coveragefn, binsize, chromSizes)
    writeWig(chromBins, chromSizes, binsize, wigfn)


main()

'''


        line = coveragefh.readline()


'''        





