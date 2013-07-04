#!/usr/bin/python
import sys
import os

def openNice(outfn, permission='U'):
    '''
    Just opens a file nicely...
    '''
    try:
        return open(outfn, permission)
    except IOError:
        print '!!!!!'
        print 'Error opening %s' % (outfn)
        print '!!!!!'
        return    

def gunzipFastqs(sampledir, readtype, species, outprefix = 'FASTQ', singleSubstr='R1', pairedSubstr='R2'):
    '''
    Unzips and concatenates the FASTQ files from an Illumina experiment
    for downstream analysis.

    Returns a dict of the names of the single and paired end file names for other functions to use.

    The concatenated file is stored in a subdirectory <outdir> of the 
    sample directory <sampledir>.

    The Substr identifies the read file as single or paired end.
    '''
    from glob import glob
    import gzip 
    outdir = sampledir + '/' + '_'.join([outprefix, readtype, species])

    print "Uncompressing files for %s to %s..." % (sampledir, outdir)
    if os.path.isdir(outdir):
        print "#####"
        print "Warning:  It looks like the unzipped FASTQ files are already present in %s" % (outdir)
        print "#####"
        soutfn = glob(outdir + '/' + '*' + singleSubstr + '.all.fastq')
        poutfn = glob(outdir + '/' + '*' + pairedSubstr + '.all.fastq')
        outfns = {}
        if soutfn:
            outfns['single'] = soutfn[0]
        else:
            return
        if poutfn:
            outfns['paired'] = poutfn[0]
        return outfns
    else:
        os.mkdir(outdir)
    
    fastqgzs = glob(sampledir + '/' + '*.gz')

    samplename = fastqgzs[0].split('/')[-1].split(singleSubstr)[0] 
    
    singlefqgzs, pairedfqgzs = [], []
    for fastqgz in fastqgzs:
        if singleSubstr in fastqgz:
            singlefqgzs.append(fastqgz)
        elif pairedSubstr in fastqgz:
            pairedfqgzs.append(pairedSubstr)
    for fastqgz in fastqgzs:
        if (fastqgz not in singlefqgzs) and (fastqgz not in pairedfqgzs):
            print "#####"
            print "WARNING:  Skipping %s because it doesn't contain a single or paired end FASTQ substring!"            
            print "#####"
    
    outfns = {}

    if len(singlefqgzs) > 0:
        soutfn = outdir + '/' + samplename + singleSubstr + '.all.fastq'
        os.system('zcat ' + ' '.join(singlefqgzs) + '> ' + soutfn)
        outfns['single'] = soutfn
    if len(pairedfqgzs) > 0:
        poutfn = outdir + '/' + samplename + pairedSubstr + '.all.fastq'
        os.system('zcat ' + ' '.join(pairedfqgzs) + '> ' + poutfn)
        outfns['paired'] = poutfn

    return outfns

def printUsage():
    '''
    Print the usage message for the program.
    '''
    print "Usage:  python ./pipeline.py <sampledir> <readtype> <species>"
    sys.exit()

def trimFastq(fastqfn, adapterseq):
    '''
    Trims given adapter sequences out of the fastq files.

    Returns the name of the trimmed fastq file, which is hard-coded to 
    end in .adaptercleaned.fastq
    '''
    import os

    fastqfh = openNice(fastqfn,'r')
    outfn = fastqfn.split('.')[0] + '.adaptercleaned.fastq'

    if os.path.isfile(outfn):
        print "#####"
        print "WARNING:  It looks like the fastq file %s has already been cleaned!" % (fastqfn)
        print "#####"
        return outfn

    outfh = openNice(outfn, 'w')
    trimcountsfn = fastqfn.split('.')[0] + '.trimcounts'
    countsfh = openNice(trimcountsfn, 'w')
    
    print "\nRemoving adapter sequences from %s..." % (fastqfn)

    # we want to count and report the number of reads and the 
    # number trimmed
    count = 0
    trimcount = 0
    while (True):
        for i in range(4):
            seqname = fastqfh.readline().strip()
            if seqname == '': 
                countsfh.write('Total Reads = %s\n' % (count))
                countsfh.write('Trimmed Reads = %s\n' % (trimcount))
                percent = float(trimcount) / count * 100
                countsfh.write('Percent trimmed = %s\n' % (percent))

                print "\tFrom a total of %s, %s were trimmed (%s percent)." % (count, trimcount, percent)

                fastqfh.close()
                outfh.close()
                countsfh.close()

                #print outfn
                return outfn

            read = fastqfh.readline().strip()
            untrimlen = len(read)
            read = removeAdapter(read, adapterseq)
            if len(read) != untrimlen:
                trimcount += 1
            # if the read is entirely composed of adapter, then leave a pair of 'Ns' to
            # preserve paired-end read order.  Two Ns are required for bowtie to handle them.
            if read == '':
                read = 'NN'
            comment = fastqfh.readline().strip()
            qscore = fastqfh.readline().strip()
            outfh.write(seqname + '\n')
            outfh.write(read + '\n')
            outfh.write(comment + '\n')
            outfh.write(qscore[0:len(read)] + '\n')
            count += 1
        if count % 1e6 == 0:
            print count, read

def removeAdapter(read, adapterseq):
    '''
    Removes an adaptersequence from a read by taking a 4-nt match and extending it.
    '''
    #print read
    idx = read.rfind(adapterseq)
    if idx != -1:
        #print 'Found the adapter at pos %s in read %s' % (idx, read)
        return read[0:idx]
    else:
        idx = read.rfind(adapterseq[0:3])
        if idx != -1:
            pad = len(read) - idx
            idx = read.rfind(adapterseq[0:pad])
            if idx != -1:
                #print 'Found the incomplete adapter at pos %s in read %s' % (idx, read)
                return read[0:idx]
            else:
                #print 'Not even incomplete adapter found in read %s' % (read)
                return read
        else:
            #print 'No adapter found in read %s' % (read)
            return read

def indexReadsBWA(fqs, readtype, species, singleSubstr='R1', pairedSubstr='R2'):
    '''
    Index the reads with BWA.
    '''
    from glob import glob
    import os
    import subprocess as SP
    from multiprocessing import cpu_count

    if species == 'mm10':
        bwaidx = '/home/matt/data/Reference/mm10/genome.fa'
    else:
        print "Error:  invalid species name for %s" % (species)
        return

    saifns = {}
    for fq in fqs:
        # make a directory to store things in
        path = os.path.dirname(fqs[fq])
        bwadir=path + '/bwa'
        if not os.path.isdir(bwadir):
            os.mkdir(bwadir)
        saifn=bwadir + '/' + os.path.basename(fqs[fq]).split('.fastq')[0] + '.sai'
        if saifn.find(singleSubstr) != -1:
            saifns['single'] = saifn
        elif saifn.find(pairedSubstr) != -1:
            saifns['paired'] = saifn
        saioutfn=bwadir + '/' + os.path.basename(fqs[fq]).split('.fastq')[0] + '.saiout'
        saioutfh = openNice(saioutfn, 'w')
        print "Making BWA idx file %s for fastq file %s" % (saifn, fqs[fq])
        args = ['bwa', 'aln', '-t', str(cpu_count()-1), bwaidx, fqs[fq], '-f', saifn] 
        proc = SP.Popen(args, stdout=saioutfh, stderr=saioutfh)
        proc.communicate()
        saioutfh.close()
    
    return saifns
        
def alignBWA(trimfqs, saifns, readtype, species, singleSubstr='R1', pairedSubstr='R2'):
    '''
    Align reads with BWA.
    '''
    import os
    from multiprocessing import cpu_count
    import subprocess as SP

    if species == 'mm10':
        bwaidx = '/home/matt/data/Reference/mm10/genome.fa'
    else:
        print "Error:  invalid species name for %s" % (species)
        return

    for saifn in saifns:
        if saifn == 'single':
            print 'Aligning single-end reads...'
            bwadir = os.path.dirname(saifns['single'])
            samfn = saifns[saifn].split('.sai')[0] + '.sam'
            samoutfn = saifns['single'].split('.sai')[0] + '.samout'
            samoutfh = openNice(samoutfn, 'w')
            print "\tAligning single-end reads with BWA to %s." % (saifns[saifn])
            print "\t\tReads from %s \n\t\twith %s for the index..." % (trimfqs['single'], saifns['single'])
            args = [ 'bwa', 'samse', bwaidx, saifns['single'], trimfqs['single'], '-f', samfn]
            proc = SP.Popen(args, stdout=samoutfh, stderr=samoutfh)
            proc.communicate()

            return samfn

        elif saifn == 'paired':
            print 'Aligning single-end reads...'
            bwadir = os.path.dirname(saifns[saifn])
            samfn = saifns['paired'].split('.sai')[0] + '.sam'
            samoutfn = saifns['paried'].split('.sai')[0] + '.samout'
            samoutfh = openNice(samoutfn, 'w')
            print "\tAligning paired-end reads with BWA to %s." % (saifns[saifn])
            print "\t\tReads from %s \n\t\twith %s for the index..." % (trimfqs['single'], saifns['single'])
            args = [ 'bwa', 'sampe', bwaidx, saifns['single'], saifns['paired'], trimfqs['single'], trimfqs['paired'], '-f', samfn]
            proc = SP.Popen(args, stdout=samoutfn, stderr=samoutfn)
            proc.communicate()

            return samfn

def makeBAM(samfn):
    '''Uses samtools to make a .bam file from a .sam file.
    '''
    # how samtools sort works...                                                                                                                                  
    import subprocess as SP
    bamfn = samfn.split('.sam')[0] + '.bam'
    print 'Converting %s to .bam format...' % (samfn)
    args = ['samtools', 'view', '-bS', samfn, '-o', bamfn]
    proc = SP.Popen(args, stdout=SP.PIPE, stderr=SP.PIPE)
    proc.wait()
    sout, serr = proc.communicate()
    print sout, serr
    print '\t...done.'

    print 'Sorting %s...' % (bamfn)
    # this split below is required b/c samtools sort wants to always add .bam to the file
    args = ['samtools', 'sort', bamfn, bamfn.split('.bam')[0]]
    proc = SP.Popen(args, stdout=SP.PIPE, stderr=SP.PIPE)
    proc.wait()
    sout, serr = proc.communicate()
    print sout, serr
    print '\t...done.'

    baifn = bamfn.split('.bam')[0] + '.bai'
    print 'Making .bai file for %s...' % (bamfn)
    args = ['samtools', 'index', bamfn, baifn]
    proc = SP.Popen(args, stdout=SP.PIPE, stderr=SP.PIPE)
    proc.wait()
    sout, serr = proc.communicate()
    print sout, serr
    print '\t...done.'

    print "Removing exact mapped duplicates..." 
    #args = ['samtools', 'rmdup', '-s', bamfn, bamfn.split('.bam')[0] + '.bam']
    rmdupfn = bamfn.split('.bam')[0] + '.rmdup'
    args = ['samtools', 'rmdup', '-s', bamfn, rmdupfn]
    proc.wait()
    proc = SP.Popen(args, stdout=SP.PIPE, stderr=SP.PIPE)
    sout, serr = proc.communicate()
    print sout, serr
    print '\t...done.'

    os.rename(bamfn, bamfn.split('.bam')[0] + '.withduplicates.bam')
    os.rename(rmdupfn, bamfn)

    print 'Making new .bai file for %s...' % (bamfn)
    args = ['samtools', 'index', bamfn, baifn]
    proc = SP.Popen(args, stdout=SP.PIPE, stderr=SP.PIPE)
    proc.wait()
    sout, serr = proc.communicate()
    print sout, serr
    print '\t...done.'

    return bamfn

def splitSAM(samfn):
    '''
    '''
    import subprocess as SP

    outfns = []
    outfns.append(samfn.split('.sam')[0] + '.watson.sam')
    outfns.append(samfn.split('.sam')[0] + '.crick.sam')
    with open(outfns[0], 'w') as watsonOutfh:
        with open(outfns[1], 'w') as crickOutfh:
            with open(samfn) as samfh:
                for line in samfh:
                    if line.startswith('@'):
                        watsonOutfh.write(line)
                        crickOutfh.write(line)
                        continue
                    # if the read is on the reverse strand, send it to crick, otherwise, watson
                    line = line.split('\t')
                    if (int(line[1]) & 16):
                        crickOutfh.write('\t'.join(line))
                    else:
                        watsonOutfh.write('\t'.join(line))

    return outfns
    
def mapBWA(trimfqs, readtype, species):
    '''
    Wrapper fucntion for mapping reads with BWA (as opposed to some other read mapper).
    '''
    import os

    testbamfn = os.path.dirname(trimfqs[readtype]) + '/' + 'bwa' + '/' + os.path.basename(trimfqs[readtype]).split('.fastq')[0] + '.bam'
    if os.path.isfile(testbamfn):
        print "######"
        print "WARNING:  BAM file already exists!  %s" % (testbamfn)
        print "######"
        return [testbamfn]

    saifns = indexReadsBWA(trimfqs, readtype, species)
    samfn = alignBWA(trimfqs, saifns, readtype, species)

    bamfns = []
    
    samfns = splitSAM(samfn)
    for samfn in samfns:
        bamfns.append(makeBAM(samfn))
    return bamfns
    bamfns.append(makeBAM(samfn))
    return bamfns

def calculateMappedLibrarySize(bamfn):
    '''
    '''
    import subprocess as SP
    
    args = ['samtools', 'idxstats', bamfn]
    proc = SP.Popen(args, stdout=SP.PIPE, stderr=SP.PIPE)
    # releveant output is from stderr
    stderr, stdout = proc.communicate()
    size = 0
    for line in stderr.split('\n'):
        line = line.split()
        if line:
            size += int(line[2])
    
    return size

def makeNormalizedBedGraph(bamfn, species, insertSize=150, scaledLibrarySize=10e7):
    '''
    '''
    import subprocess as SP
    import gzip

    #if bamfn.endswith('.mm10.bam'):
    chromsizefn = 'mm10.chrom.sizes'
    species = 'mm10'
    #elif bamfn.endswith('.dpse.bam'):
    #    chromsizefn = 'dp3.chrom.sizes'
    #    species = 'dp3'
    #elif species == 'dm3':
    #    chromsizefn = 'dm3.chrom.sizes'
    #elif species == 'dp3':
    #    chromsizefn = 'dp3.chrom.sizes'
    #else:
    #    print "Error:  invalid species name for %s" % (species)
    #    return
    chromsizes = {}
    with openNice(chromsizefn, 'r') as chromsizefh:
        for line in chromsizefh:
            line = line.strip().split('\t')
            if line == ['']: 
                continue
            else:
                print line
                chromsizes[line[0]] = int(line[1])

    print "Making Normalized BedGraph for %s..." % (bamfn)
    
    librarySize = calculateMappedLibrarySize(bamfn)
    print "There are %s reads in the mapped library." % (librarySize)
    with openNice(bamfn.split('.bam')[0] + '.mappedLibrarySize', 'w') as sizefh:
        sizefh.write(str(librarySize))

    # write bed file from bam file
    bedfn = bamfn.split('.bam')[0] + '.bed'
    bedfh = openNice(bedfn, 'w')
    args = ['bamToBed', '-i', bamfn]
    proc = SP.Popen(args, stdout=SP.PIPE, stderr=SP.PIPE)
    bedPIPE = proc.communicate()[0].split('\n')
    
    count = 0

    for line in bedPIPE:
        count += 1
        line = line.split()
        # somewhere the bamToBed program introduces a blank line, this handles it
        try:
            chrom = line[0]
        except IndexError as e:
            #print line, count, e
            continue
        start = int(line[1])
        stop = int(line[2])
        strand = line[5]
        readlen = stop - start
        if strand == '+':
            stop = start + insertSize
            if stop > chromsizes[chrom]:
                stop = chromsizes[chrom]
        elif strand == '-':
            stop = stop + readlen
            start = stop - insertSize
            if start < 1:
                start = 1
            if stop > chromsizes[chrom]:
                stop = chromsizes[chrom]
        bedfh.write('\t'.join([chrom,str(start),str(stop),strand]) + '\n')
    bedfh.close()

    # use GNU sort to sort the bed file
    print "Sorting the bed file %s with GNU sort..." % (bedfn)
    args = ['sort', '-k1,1', '-k2,2n', bedfn, '-o', bedfn]
    proc = SP.Popen(args, stdout=SP.PIPE, stderr=SP.PIPE)
    sout, serr = proc.communicate()
    if serr: print serr
    
    # make the density coverage file
    wigfn = bamfn.split('.bam')[0] + '.wig'
    wigfh = open(wigfn, 'w')

    print "Making the density coverage file with genomeCoverageBed..."
    args = ['genomeCoverageBed', '-i', bedfn, '-g', chromsizefn, '-dz']
    proc = SP.Popen(args, stdout=wigfh, stderr=wigfh)
    proc.wait()
    #densities = proc.communicate()[0].split('\n')
    wigfh.close()
    
    #print "density file written, exiting program..."

    print "Writing scaled density file..."
    wigfh = open(wigfn, 'r')
    scaledwigfn = bamfn.split('.bam')[0] + '.scaled.wig'
    scaledwigfh = open(scaledwigfn, 'w')
    chrom = ''
    lastline = ''
    for line in wigfh:
        line = line.split('\t')
        try:
            curchrom = line[0]
            # wig files are 1-based
            pos = int(line[1]) + 1
            count = float(line[2]) * scaledLibrarySize / librarySize
        except Exception as e:
            continue
        if not curchrom == chrom:
            print "\t...%s..." % (curchrom)
            chrom = curchrom
            scaledwigfh.write('variableStep chrom=%s\n' % (chrom))
        outline = '\t'.join([str(pos), str(count)])
        scaledwigfh.write(outline + '\n')
        lastline = line
    scaledwigfh.close()
    wigfh.close()

    # Create BigWig 
    bwfn = bamfn.split('.bam')[0] + '.bw'
    print "Creating BigWig file %s..." % (bwfn)
    args = ['wigToBigWig', scaledwigfn, chromsizefn, bwfn]
    proc = SP.Popen(args)
    proc.wait()

    # Create BedGraph
    bgfn = bamfn.split('.bam')[0] + '.bedGraph'
    print "Creating BedGraph file %s..." % (bgfn)
    args = ['bigWigToBedGraph', bwfn, bgfn]
    print args
    proc = SP.Popen(args)
    proc.wait()

def cleanupSeqFiles(fqs, bamfn, keepsam = False):
    '''
    '''
    import os
    from glob import glob

    prefix = os.path.dirname(fqs['single'])
    fns = glob(prefix + '/*.fastq')
    for fn in fns: 
        os.remove(fn)

    prefix += '/bwa/'
    fns = []
    fns.extend(glob(prefix + '/*.sai'))
    fns.extend(glob(prefix + '/*.bw'))
    fns.extend(glob(prefix + '/*.bed'))
    fns.extend(glob(prefix + '/*.wig'))

    for fn in fns:
        os.remove(fn)

    samfns = glob(prefix + '/*.sam')

    if keepsam == True:
        for samfn in samfns:
            os.system("gzip %s" % (samfn))
    else:
        for samfn in samfns:
            os.remove(samfn)


if __name__ == '__main__':

    if len(sys.argv) != 4:
        printUsage()
        
    sampledir = sys.argv[1]
    readtype = sys.argv[2]
    species = sys.argv[3]
    adapterSE='GATCGGAAGAGCACACGTCTGAACTCCAGTCAC'
    adapterPE='ATAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT'

    # gunzip and concatenate the fq.gz files
    # fqs is a dictionary with the concatenated filename for 
    # single and paired files
    fqs = gunzipFastqs(sampledir, readtype, species)
    if not fqs:
        print "Exiting..."
        sys.exit()

    # trim the adapter sequences from the reads
    # capture the filenames in the return to pass
    # to the mapping function of your choice
    trimfqs = {}
    for fq in fqs:
        if fq == 'single':
            trimfqs['single'] = trimFastq(fqs[fq], adapterSE)
        elif fq == 'paried':
            trimfqs['paired'] = trimFastq(fqs[fq], adapterPE)
            
    bamfns = mapBWA(trimfqs, readtype, species)
    for bamfn in bamfns:
        makeNormalizedBedGraph(bamfn, species, insertSize=350)
        cleanupSeqFiles(fqs, bamfn, keepsam=False)
