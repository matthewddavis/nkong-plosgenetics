#!/g/software/bin/python-2.7
import sys
import logging
import os
from readTools import gunzipFastqs
from littleTools import openNice

def openLogger(log_fn, log_name):
    # clear any previous log file of the same name
    if os.path.isfile(log_fn):
        os.remove(log_fn)
    # log for pipeline
    logger = logging.getLogger(log_name)
    logger.setLevel(logging.DEBUG)
    # log file captures DEBUG and above
    fh = logging.FileHandler(log_fn)
    fh.setLevel(logging.DEBUG)
    # console log is ERROR and above
    ch = logging.StreamHandler()
    ch.setLevel(logging.ERROR)
    # create formatter and add it to the handlers
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)
    # add the handlers to the logger
    logger.addHandler(fh)
    logger.addHandler(ch)

    return logger

def trimFastq(fastq_fn, adapter_seq):
    '''
    Trims given adapter sequences out of the fastq files.

    Returns the name of the trimmed fastq file, which is hard-coded to 
    end in .adaptercleaned.fastq
    '''
    import os

    fastqfh = openNice(fastq_fn,'r')
    outfn = fastq_fn.split('.')[0] + '.adaptercleaned.fastq'

    if os.path.isfile(out_fn):
        print "#####"
        print "WARNING:  It looks like the fastq file %s has already been cleaned!" % (fastq_fn)
        print "#####"
        return out_fn

    out_fh = openNice(out_fn, 'w')
    trimcounts_fn = fastq_fn.split('.')[0] + '.trimcounts'
    counts_fh = openNice(trimcounts_fn, 'w')
    
    print "\nRemoving adapter sequences from %s..." % (fastq_fn)

    # we want to count and report the number of reads and the 
    # number trimmed
    count = 0
    trim_count = 0
    while (True):
        for i in range(4):
            seq_name = fastq_fh.readline().strip()
            if seq_name == '': 
                counts_fh.write('Total Reads = %s\n' % (count))
                counts_fh.write('Trimmed Reads = %s\n' % (trim_count))
                percent = float(trim_count) / count * 100
                counts_fh.write('Percent trimmed = %s\n' % (percent))

                print "\tFrom a total of %s, %s were trimmed (%s percent)." % (count, trim_count, percent)

                fastq_fh.close()
                out_fh.close()
                counts_fh.close()

                return out_fn

            read = fastq_fh.readline().strip()
            untrim_len = len(read)
            read = removeAdapter(read, adapter_seq)
            if len(read) != untrim_len:
                trim_count += 1
            # if the read is entirely composed of adapter, then leave a pair of 'Ns' to
            # preserve paired-end read order.  Two Ns are required for bowtie to handle them.
            if read == '':
                read = 'NN'
            comment = fastq_fh.readline().strip()
            qscore = fastq_fh.readline().strip()
            out_fh.write(seq_name + '\n')
            out_fh.write(read + '\n')
            out_fh.write(comment + '\n')
            out_fh.write(qscore[0:len(read)] + '\n')
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

def makeNormalizedBedGraph(bamfn, species='mm10', insertSize=150, scaledLibrarySize=10e7):
    '''
    '''
    import subprocess as SP
    import gzip

    chromsizefn = 'reference/mm10.chrom.sizes'

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

def runMacs14(treat_fn, species='mm10', control_fn='', macs_bin='/g/software/bin/macs14'):
    '''                                                                                                                                                                
    - treatfn should be a processed and sorted bamfn                                                                                                                   
    - if controlfn is unspecified, MACS14 will use the no input control model                                                                                          
    '''
    if species == 'mm10':
        genome_size = '1.87e9'
    else:
        print "Invalid species type %s" % (species)
        raise Error

    data_name = treat_fn.split('.bam')[0] + '.macs'

    if control_fn:
        print "\tRunning MACS14 with %s as the treatment and %s as the control..." % (treat_fn, control_fn)
        args = [macs_bin, '-t', treat_fn, '-c', control_fn, '-n', data_name, '-fBAM', '-g', genome_size]
    else:
        print "\tRunning MACS14 with %s as the treatment and no control file..." % (treat_fn)
        args = [macs_bin, '-t', treat_fn, '-n', data_name, '-fBAM', '-g', genome_size]

    proc = SP.Popen(args, stdout=SP.PIPE, stderr=SP.PIPE)
    proc.wait()
    sout, serr = proc.communicate()
    print sout, serr
    print '\t...done.'

class SeqSample(object):
    '''
    '''
    def __init__(self, sample_dir, log_name='pipeline'):
        self.logger = logging.getLogger('pipeline')
        
        self.attrs = {}
        self.attrs['sample_dir'] = sample_dir
        self.loadSampleAttrs()

    def loadSampleAttrs(self):
        '''
        Each sample has a sample_info file in the sample_dir.
        Each line of this file is a semi-colon separated attr
        key:value pair with whitespace allowed.

        This function parses these attrs into the SeqSample object.
        '''
        sample_info_fn = self.attrs['sample_dir'] + '/' + 'sample_info'
        for line in open(sample_info_fn):
            if line.startswith('#'):
                continue
            line = line.split(';')
            line = [i.strip() for i in line]
            self.attrs[line[0]] = line[1]

        logger.info('Sample loaded:')
        for attr in self.attrs:
            logger.info('\t%s : %s' % (attr, self.attrs[attr]))



if __name__ == '__main__':


    logger = openLogger(log_fn='pipeline.log', log_name='pipeline')

    sample = SeqSample(sample_dir='data/Sample_NKRT1/')

    sys.exit()


    sample_dir = sys.argv[1]
    adapterSE='GATCGGAAGAGCACACGTCTGAACTCCAGTCAC'

    # gunzip and concatenate the fq.gz files
    # fqs is a dictionary with the concatenated filename for 
    # single and paired files
    fqs = gunzipFastqs(sample_dir, read_type='single', species='mm10')
    if fqs:
        # trim the adapter sequences from the reads
        # capture the filenames in the return to pass
        # to the mapping function of your choice
        trim_fqs = {}
        for fq in fqs:
            if fq == 'single':
                trimfqs['single'] = trimFastq(fqs[fq], adapterSE)
            elif fq == 'paried':
                trimfqs['paired'] = trimFastq(fqs[fq], adapterPE)
        fqs = trimfqs
    
    bamfns = mapBWA(fqs, readtype, species)
    for bamfn in bamfns:
        makeNormalizedBedGraph(bamfn, species, insertSize=350)
        cleanupSeqFiles(fqs, bamfn, keepsam=False)

        runMacs14(treat_fn=bamfn, species='mm10', control_fn='', macs_bin='/g/software/bin/macs14')
