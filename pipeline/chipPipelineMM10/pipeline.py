#!/g/software/bin/python-2.7
import sys
import logging
import os
from readTools import gunzipFastqs
from littleTools import openNice
import inspect

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
    ch.setLevel(logging.DEBUG)
    # ch.setLevel(logging.INFO)
    # create formatter and add it to the handlers
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)
    # add the handlers to the logger
    logger.addHandler(fh)
    logger.addHandler(ch)

    return logger

def parseOptions():
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("--topDataDir", dest="top_data_dir",
                      help="this is the top-level dir containing subdirs of sequencing reads", metavar="DIR")

    (options, args) = parser.parse_args()
    return (options, args)






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


    
class SeqSample(object):
    '''
    '''
    def __init__(self, sample_dir, log_name='pipeline', cpu_count=1, single_str='R1', paired_str='R2'):
        self.logger = logging.getLogger('pipeline')

        self.sample_dir = sample_dir
        self.cpu_count = str(cpu_count)
        self.single_str = single_str
        self.paired_str = paired_str
        
        self.attrs = {}
        self.loadSampleAttrs()

        self.use_untrimmed = None

    def loadSampleAttrs(self):
        '''
        Each sample has a sample_info file in the sample_dir.
        Each line of this file is a semi-colon separated attr
        key:value pair with whitespace allowed in the values.

        E.g.
        sample_id; Sample_NKRT1
        common_name; Mef2 ChIP-exo
        data_type; chip-exo
        read_type; single
        read_length; 50
        adapter_seq; GATCGGAAGAGCACACGTCTGAACTCCAGTCAC

        This function parses these attrs into the SeqSample object.
        '''
        sample_info_fn = self.sample_dir + '/' + 'sample_info'
        for line in open(sample_info_fn):
            if line.startswith('#'):
                continue
            line = line.split(';')
            line = [i.strip() for i in line]
            self.attrs[line[0]] = line[1]

        self.logger.info('Sample loaded:')
        for attr in self.attrs:
            self.logger.info('\t%s : %s' % (attr, self.attrs[attr]))



    def gunzipFastqs(self):
        '''
        Unzips and concatenates the FASTQ files from an Illumina experiment
        for downstream analysis.

        The concatenated file is stored in a subdirectory <outdir> of the 
        sample directory <sampledir>.

        The single_str or paired_str identifies the read file as single or paired end.
        For VCGSL reads, this is R1 for single and R2 for paired.
        '''
        from glob import glob
        import gzip
        from datetime import datetime
        from multiprocessing import Pool

        def _checkForFastq(out_dir):
            '''
            Check to see if the uncompressed files already exist from a previous run.

            The out_dir can be anywhere, but usually is the sample_dir.
            '''
            read_type = self.attrs['read_type']
            out_fns = glob(out_dir + '/' + '*' + '.all.fastq')
            if len(out_fns) == 0:
                return False

            # There is at least one all.fastq, but there should be at most two all.fastq files
            # (in the case of paired-end reads)
            assert len(out_fns) <= 2
            for out_fn in out_fns:
                if read_type == 'single':
                    # check for single-end reads of any length
                    if self.single_str in out_fn:
                        self.logger.warn("An uncompressed single-end fastq file for %s already exists: %s" % (self.attrs['sample_id'], out_fn) )
                        self.logger.warn("\t...skipping decompression.")
                        self.single_fastq_fn = out_fn
                        self.paired_fastq_fn = None
                        return True
                elif read_type == 'paired':
                    # these are paired-end reads of any length
                    if self.single_str in out_fn:
                        self.single_fastq_fn = out_fn
                    elif self.paired_str in out_fn:
                        self.paired_fastq_fn = out_fn
            # make sure both single-end and paired-end were
            # found, if this is a paired-end sample
            if self.single_fastq_fn and self.paired_fastq_fn:
                self.logger.warn("Uncompressed paired-end fastq files for %s already exist: %s" % (self.attrs['sample_id'], out_fn) )
                self.logger.warn("\t...skipping decompression.")
                return True
            else:
                return False

        def _sortFastqgzFiles():
            '''
            Separates the reads into single and paired, if need be.
            '''
            fastq_gzs = glob(self.sample_dir + '/' + '*.gz')
            single_fq_gzs, paired_fq_gzs = [], []
            for fastq_gz in fastq_gzs:
                if self.single_str in fastq_gz:
                    single_fq_gzs.append(fastq_gz)
                elif self.paired_str in fastq_gz:
                    paired_fq_gzs.append(fastq_gz)
                else:
                    self.logger.warn("Skipping %s because it doesn't contain a single-end or paired-end FASTQ substring." % (fastq_gz) )

            if self.attrs['read_type'] == 'single':
                assert len(single_fq_gzs) > 0
            elif self.attrs['read_type'] == 'paired':
                assert len(single_fq_gzs) > 0
                assert len(paired_fq_gzs) > 0
            else:
                self.logger.error('!!! The read_type is not present or not interpretable in the sample_info file.')
                self.logger.error('!!! Exiting on line %s...' % (inspect.currentframe().f_lineno) )
                sys.exit()
            return single_fq_gzs, paired_fq_gzs

        ### Begin top-level gunzipFastqs()

        self.single_fastq_fn = None
        self.paired_fastq_fn = None
        
        # if the fastq file(s) exist already, return
        if _checkForFastq(self.sample_dir):
            return

        single_fq_gzs, paired_fq_gzs = _sortFastqgzFiles()

        if self.attrs['read_type'] == 'single':
            self.single_fastq_fn = self.sample_dir + '/' + self.attrs['sample_id'] + '_' + self.single_str + '.all.fastq'
            self.paired_fastq_fn = None
            self.logger.info("Uncompressing single-end fastq files for %s..." % (self.sample_dir))
            cmd = 'zcat ' + ' '.join(single_fq_gzs) + '> ' + self.single_fastq_fn
            os.system(cmd)
        if self.attrs['read_type'] == 'paired':
            self.single_fastq_fn = self.sample_dir + '/' + self.attrs['sample_id'] + '_' + self.single_str + '.single.all.fastq'
            self.paired_fastq_fn = self.sample_dir + '/' + self.attrs['sample_id'] + '_' + self.paired_str + '.paired.all.fastq'
            self.logger.info("Uncompressing single-end fastq files for %s..." % (self.sample_dir))
            cmd = 'zcat ' + ' '.join(single_fq_gzs) + '> ' + self.single_fastq_fn
            os.system(cmd)
            self.logger.info("Uncompressing paired-end fastq files for %s..." % (self.sample_dir))
            cmd = 'zcat ' + ' '.join(paired_fq_gzs) + '> ' + self.paired_fastq_fn
            os.system(cmd)

    def trimFastqFiles(self, keep_untrimmed_fastq=False):
        '''
        Trims given adapter sequences out of the uncompressed fastq file.

        Unless keep_untrimmed_fastq is set to True, the *.all.fastq files
        are deleted to conserve disk space.
        '''
        def _checkForTrimmedFastq():
            '''
            Checks to see if the trimmed file(s) exist from a previous run.
            '''
            if self.attrs['read_type'] == 'single':
                if os.path.isfile(single_trimmed_fastq_fn):
                    self.logger.warn("It looks like the fastq file %s has already been trimmed!" % (self.single_fastq_fn) )
                    return True
            if self.attrs['read_type'] == 'paired':
                if os.path.isfile(single_trimmed_fastq_fn) and os.path.isfile(paired_trimmed_fastq_fn):
                    self.logger.warn("It looks like the fastq file %s has already been trimmed!" % (self.single_fastq_fn) )
                    self.single_trimmed_fastq_fn = single_trimmed_fastq_fn
                    self.logger.warn("It looks like the fastq file %s has already been trimmed!" % (self.paired_fastq_fn) )
                    self.paired_trimmed_fastq_fn = paired_trimmed_fastq_fn
                    return True
                
        def _removeUntrimmedFastqs():
            if self.single_fastq_fn != None:
                os.remove(self.single_fastq_fn)
            if self.paired_fastq_fn != None:
                os.remove(self.paired_fastq_fn)

        def _trimFastqFile(fastq_fn):
            '''
            Trims the full or partial adapter sequences from the reads
            in the specified fastq file.
            '''
            def _removeAdapter(read, adapter_seq):
                '''
                Removes an adapter sequence from a read by taking a 4-nt match and extending it.
                '''
                idx = read.rfind(adapter_seq)
                if idx != -1:
                    self.logged.debug('Found the adapter at pos %s in read %s' % (idx, read))
                    return read[0:idx]
                else:
                    idx = read.rfind(adapter_seq[0:3])
                    if idx != -1:
                        pad = len(read) - idx
                        idx = read.rfind(adapter_seq[0:pad])
                        if idx != -1:
                            self.logger.debug('Found the incomplete adapter at pos %s in read %s' % (idx, read) )
                            return read[0:idx]
                        else:
                            self.logger.debug('Not even incomplete adapter found in read %s' % (read) )
                            return read
                    else:
                        self.logger.debug('No adapter found in read %s' % (read) )
                return read

            ### _trimFastqFile body starts here
            fastq_fh = open(fastq_fn, 'r')
            
            trimmed_fn = fastq_fn.replace('.all.', '.trimmed.')
            trimmed_fh = open(trimmed_fn, 'w')
            
            trim_counts_fn = fastq_fn.split('.')[0] + '.trimcounts'
            trim_counts_fh = open(trim_counts_fn, 'w')

            self.logger.info("Removing adapter sequences from %s..." % (fastq_fn) )
            # we want to count and report the number of reads and the 
            # number trimmed
            count = 0
            trim_count = 0
            while (True):
                for i in range(4):
                    seq_name = fastq_fh.readline().strip()
                    # if the read is blank, we clean up and exit
                    if seq_name == '': 
                        trim_counts_fh.write('Total Reads = %s\n' % (count))
                        trim_counts_fh.write('Trimmed Reads = %s\n' % (trim_count))
                        percent = float(trim_count) / count * 100
                        trim_counts_fh.write('Percent trimmed = %s\n' % (percent))
                        self.logger.info("\tFrom a total of %s, %s were trimmed (%s percent)." % (count, trim_count, percent) )

                        if self.single_str in trimmed_fn:
                            self.single_trimmed_fastq_fn = trimmed_fn
                            return
                        elif self.paired_str in trimmed_fn:
                            self.paired_trimmed_fastq_fn = trimmed_fn
                            return
                        else:
                            self.logger.error('!!! The read_type is not present or not interpretable in the sample_info file.')
                            self.logger.error('!!! Exiting on line %s...' % (inspect.currentframe().f_lineno) )
                            sys.exit()

                    read = fastq_fh.readline().strip()
                    untrim_len = len(read)
                    read = _removeAdapter(read, self.attrs['adapter_seq'])
                    if len(read) != untrim_len:
                        trim_count += 1
                    # if the read is entirely composed of adapter, then leave a pair of 'Ns' to
                    # preserve paired-end read order.  Two Ns are required for bowtie to handle them,
                    # and even though the pipeline uses BWA, these files could be easily used for
                    # other read mappers
                    if read == '':
                        read = 'NN'
                    comment = fastq_fh.readline().strip()
                    qscore = fastq_fh.readline().strip()
                    trimmed_fh.write(seq_name + '\n')
                    trimmed_fh.write(read + '\n')
                    trimmed_fh.write(comment + '\n')
                    trimmed_fh.write(qscore[0:len(read)] + '\n')
                    count += 1
                if count % 1e5 == 0:
                    self.logger.info("Reads counted: %s" % (count) )
                    self.logger.info("Last read: %s" % (read) )

        ## trimFastqFiles() body starts here
        if self.attrs['read_type'] == 'single':
            single_fastq_fh = open(self.single_fastq_fn, 'r')
            single_trimmed_fastq_fn = self.single_fastq_fn.replace('.all.', '.trimmed.')
            if _checkForTrimmedFastq():
                self.single_trimmed_fastq_fn = single_trimmed_fastq_fn
                self.paired_trimmed_fastq_fn = None
                return
            else:
                _trimFastqFile(self.single_fastq_fn)

        elif self.attrs['read_type'] == 'paired':
            single_fastq_fh = open(self.single_fastq_fn, 'r')
            single_trimmed_fastq_fn = self.single_fastq_fn.replace('.all.', '.trimmed.')
            paired_fastq_fh = open(self.paired_fastq_fn, 'r')
            paired_trimmed_fastq_fn = self.paired_fastq_fn.replace('.all.', '.trimmed.')            
            if _checkForTrimmedFastq():
                self.single_trimmed_fastq_fn = single_trimmed_fastq_fn
                self.paired_trimmed_fastq_fn = paired_trimmed_fastq_fn
                return
            else:
                _trimFastqFile(self.single_fastq_fn)
                _trimFastqFile(self.paired_fastq_fn)
        
        else:
            self.logger.error('!!! The read_type is not present or not interpretable in the sample_info file.')
            self.logger.error('!!! Exiting on line %s...' % (inspect.currentframe().f_lineno) )
            import pdb
            pdb.set_trace()
            sys.exit()

        if keep_untrimmed_fastq == False:
            _removeUntrimmedFastqs()

    def indexReadsBWA(self):
        '''
        Index the reads with BWA.
        '''
        import os
        import subprocess as SP

        def _createBWADir():
            '''
            Create bwa_dir, if it doesn't already exist.
            '''
            bwa_dir = self.sample_dir + '/bwa'
            self.bwa_dir = bwa_dir

            if os.path.isdir(bwa_dir):
                self.logger.info('BWA directory already exists for sample %s.' % (self.attrs['sample_id']) )
            else:
                self.logger.info('Creating BWA directory for sample %s.' % (self.attrs['sample_id']) )
                os.mkdir(bwa_dir)

        _createBWADir()
        bwa_dir = self.bwa_dir

        if self.use_untrimmed:
            fq_fns = [self.single_fastq_fn, self.paired_fastq_fn]
        else:
            fq_fns = [self.single_trimmed_fastq_fn, self.paired_trimmed_fastq_fn]
        
        for fq_fn in [self.single_fastq_fn, self.paired_fastq_fn]:
            if fq_fn == None:
                continue
            else:
                self.logger.info("Making BWA index file for fastq file %s (sample %s)..." % (fq_fn, self.attrs['sample_id']) )
            # make a directory to store things in
            sai_fn = bwa_dir + '/' + fq_fn.replace(self.sample_dir + '/', '').replace('.fastq', '.sai')
            # if the file is already there, skip to the next file
            if os.path.isfile(sai_fn):
                continue

            if self.single_str in sai_fn:
                self.single_sai_fn = sai_fn
            elif self.paired_str in sai_fn:
                self.paired_sai_fn = sai_fn
            else:
                self.logger.error('!!! The read_type is not present or not interpretable in the sample_info file.')
                self.logger.error('!!! Exiting on line %s...' % (inspect.currentframe().f_lineno) )

            sai_log_fn = sai_fn.replace('.fastq', '.sai.log')
            sai_log_fh = open(sai_log_fn, 'w')
            args = ['bwa', 'aln', '-t', str(self.cpu_count), self.bwa_idx, fq_fn, '-f', sai_fn] 
            proc = SP.Popen(args, stdout=sai_log_fh, stderr=sai_log_fh)
            proc.communicate()
            sai_log_fh.close()
            self.logger.info("Finished making BWA idx file %s for fastq file %s (%s)." % (sai_fn, fq_fn, self.attrs['sample_id']) )

    def alignBWA(self,):
        '''
        Align reads with BWA.
        '''
        import os
        import subprocess as SP

        bwa_idx = self.bwa_idx
        bwa_dir = self.bwa_dir

        if self.use_untrimmed:
            fqs = [self.single_fastq_fn, self.paired_fastq_fn]
        else:
            fqs = [self.single_trimmed_fastq_fn, self.paired_trimmed_fastq_fn]
        
        if self.attrs['read_type'] == 'single':
            self.logger.info('Aligning single-end reads for sample %s' % (self.attrs['sample_id']) )
            sam_fn = self.single_sai_fn.replace('.sai', '.sam')
            sam_log_fn = self.single_sai_fn.replace('.sai', '.sam.log')
            sam_log_fh = open(sam_log_fn, 'w')
            self.logger.info("Aligning single-end reads for sample %s with BWA to %s." % (self.attrs['sample_id'], self.single_sai_fn) )
            args = [ 'bwa', 'samse', bwa_idx, self.single_sai_fn, fqs[0], '-f', sam_fn]
            proc = SP.Popen(args, stdout=sam_log_fh, stderr=sam_log_fh)
            proc.communicate()

        elif self.attrs['read_type'] == 'paired':
            # make sure we have a paired fastq file
            assert fqs[1] != None
            self.logger.info('Aligning paired-end reads for sample %s' % (self.attrs['sample_id']) )
            # doesn't matter which sai file you pull the name from,
            # as long as the both exist (and this implicitly checks for them)
            sam_fn = self.paired_sai_fn.replace('.sai', '.sam')
            sam_log_fn = self.paired_sai_fn.replace('.sai', '.sam.log')
            sam_log_fh = open(sam_log_fn, 'w')
            self.logger.info("Aligning paired-end reads for sample %s with BWA to %s." % (self.attrs['sample_id'], self.paired_sai_fn) )
            args = [ 'bwa', 'sampe', bwa_idx, self.single_sai_fn, self.paired_sai_fn, fqs[0], fqs[1], '-f', sam_fn]
            proc = SP.Popen(args, stdout=sam_log_fn, stderr=sam_log_fn)
            proc.communicate()

        self.sam_fn = sam_fn

    def splitSAMByStrand(self):
        '''
        For data types such as ChIP-Exo, we may want to split the reads into two sam files,
        one for each strand. This method does this using the bit flags in the sam files.

        This is typically just done for visualization in a genome browser, but could also
        be used for a custom analysis of strand biased enrichment.
        '''
        import subprocess as SP

        split_sam_watson_fn = self.sam_fn.replace('.sam', '.watson.sam')
        split_sam_crick_fn = self.sam_fn.replace('sam', 'crick.sam')
        
        split_sam_watson_fh = open(split_sam_watson_fn, 'w')
        split_sam_crick_fh = open(split_sam_crick_fn, 'w')
        sam_fh = open(self.sam_fn, 'r')

        self.logger.info("Splitting the sam file for sample %s into separate strands." % (self.attrs['sample_id']) )
        for line in sam_fh:
            # the sam header lines start with @
            if line.startswith('@'):
                split_sam_watson_fh.write(line)
                split_sam_crick_fh.write(line)
                continue
            # if the read is on the reverse strand, demarcated by the bitwise 16 operation,
            # send it to crick, otherwise, send it to watson
            line = line.split('\t')
            if (int(line[1]) & 16):
                split_sam_crick_fn.write('\t'.join(line))
            else:
                split_sam_watson_fn.write('\t'.join(line))

        self.split_sam_watson_fn = split_sam_watson_fn
        self.split_sam_crick_fn = split_sam_crick_fn

    def mapBWA(self, use_untrimmed=False):
        '''
        Mapping reads with BWA (as opposed to some other read mapper).
        '''
        import os

        def _setUntrimmedFlag():
            '''
            If use_untrimmed is set True, the untrimmed FASTQ files will be mapped.
            '''
            if use_untrimmed == True:
                assert os.path.isfile(self.single_fastq_fn)
            self.use_untrimmed = use_untrimmed
            
        def _getBWAIDX_fn(species):
            if self.attrs['species'] == 'mm10':
                bwa_idx = 'reference/mm10/genome.fa'
                self.bwa_idx = bwa_idx
            else:
                self.logger.error("Error: invalid species name for %s" % (self.attrs['species']) )
                self.logger.error('!!! Exiting on line %s...' % (inspect.currentframe().f_lineno) )
                sys.exit()

        _getBWAIDX_fn(self.attrs['species'])
        bwa_idx = self.bwa_idx
        
        # map and align the reads
        self.indexReadsBWA()
        self.alignBWA()

        # make the BAM files
        if self.attrs['data_type'] == 'chip-exo':
            splitSAMByStrand()
            self.bam_fn = samToBam(self.sam_fn)
            self.split_bam_waston_fn = samToBam(self.split_sam_watson_fn)
            self.split_bam_crick_fn =  samToBam(self.split_sam_crick_fn)
        elif self.attrs['data_type'] == 'chip-seq':
            if self.remove_duplicate_reads == True:
                self.duplicates_removed_bam_fn = removeDuplicateReads(self.bam_fn)
                    
    def samToBam(sam_fn):
        '''
        Converts sam to bam and returns the bam_fn so that any sam_fn can
        be converted and stored externally.
        '''
        import subprocess as SP

        def _sortBAM(bam_fn):
            '''
            Sorts a bam file.
            '''
            self.logger.info('Sorting bam file %s...' % (bam_fn) )
            # the replace() below is required b/c samtools sort wants to always add .bam to the file
            args = ['samtools', 'sort', bam_fn, bam_fn.replace('.bam', '')]
            proc = SP.Popen(args, stdout=SP.PIPE, stderr=SP.PIPE)
            proc.wait()
            sout, serr = proc.communicate()
            self.logger.debug(sout)
            self.logger.debug(serr)
            self.logger.info('Finished sorting bam file %s.' % (bam_fn) )

        ## SamToBam main body begins here
        bam_fn = sam_fn.replace('.sam', '.bam')
        self.logger.info('Converting %s to .bam format...' % (sam_fn) )
        args = ['samtools', 'view', '-bS', samfn, '-o', bamfn]
        proc = SP.Popen(args, stdout=SP.PIPE, stderr=SP.PIPE)
        proc.wait()
        sout, serr = proc.communicate()
        self.logger.debug(sout)
        self.logger.debug(serr)
        self.logger.info('Finished converting %s to .bam format...' % (sam_fn) )

        _sortBAM(bam_fn)
        self.makeBAI(self, bam_fn)
        
        return bam_fn

    def makeBAI(self, bam_fn):
        '''
        Creates the bam idx file required to read the .bam file.
        '''
        bai_fn = bam_fn.replace('.bam', '.bai')
        self.logger.info('Making .bai file for %s...' % (sample_) )
        args = ['samtools', 'index', bam_fn, bai_fn]
        proc = SP.Popen(args, stdout=SP.PIPE, stderr=SP.PIPE)
        proc.wait()
        sout, serr = proc.communicate()
        self.logger.debug(sout)
        self.logger.debug(serr)
        self.logger.info('Finish making .bai file for %s.' % (bam_fn) )
            
    def removeDuplicateReads(bam_fn):
        '''
        In some cases, such as ChIP-Seq, we may want to remove exact
        duplicate reads from the samples.
        '''
        self.logger.info("Removing exact mapped duplicates for sample %s from bam file %s..." % (self.attrs['sample_id'], bam_fn) )
        duplicates_removed_bam_fn = bamfn.replace('.bam', '.rmdup.bam')
        args = ['samtools', 'rmdup', '-s', bamfn, duplicates_removed_bam_fn]
        proc = SP.Popen(args, stdout=SP.PIPE, stderr=SP.PIPE)
        proc.wait()
        sout, serr = proc.communicate()
        self.logger.debug(sout)
        self.logger.debug(serr)
        self.logger.info('Finished removing exact mapped duplicates for sample %s from bam file %s' % (self.attrs['sample_id'], bam_fn) )

        # need to remake the .bai file after removing duplicate reads
        self.logger.info('Re-making the .bai file for %s after removing duplicate reads' % (self.attrs['sample_id']) )
        self.makeBAI(duplicates_removed_bam_fn)
        
        return duplicates_removed_bam_fn

    def runMACS(self, treatment_bam_fn, control_sample_path=None, macs_bin='/usr/bin/macs14'):
        '''
        Uses MACS to call peaks in a dataset.

        The treatment_bam_fn can must be a sorted bam file.

        If control_sample_path is unspecified, MACS14 will use the no input control model                                                                                          
        '''
        import os

        def _getControlBamFn():
            '''
            Gets the control bam filename for the control sample path specified in the
            sample_info file.
            '''
            from glob import glob
            bam_fns = glob(control_sample_path + '/' + self.bwa_dir + '/' + '*.bam')
            # need to add some code here for the split read bam files (15 July 2013)
            assert len(bam_fns) == 1
            return bam_fns[0]
        

        # make sure MACS is installed
        assert os.path.isfile(macs_bin)

        if self.attrs['species'] == 'mm10':
            genome_size = '1.87e9'
        else:
            self.logger.error('!!! The species is not specified or not interpretable in the sample_info file.')
            self.logger.error('!!! Exiting on line %s...' % (inspect.currentframe().f_lineno) )
            sys.exit()

        macs_fn = treatment_bam_fn.replace('.bam', '.macs')

        if control_sample_path == None:
            self.logger.info("Running MACS with %s as the treatment and no control file..." % (treatment_bam_fn) )
            args = [macs_bin, '-t', treatment_bam_fn, '-n', macs_fn, '-fBAM', '-g', genome_size]
        else:
            control_bam_fn = _getControlBamFn()
            print "\tRunning MACS14 with %s as the treatment and %s as the control..." % (treatment_bam_fn, control_bam_fn)
            args = [macs_bin, '-t', treatment_bam_fn, '-c', control_bam_fn, '-n', macs_fn, '-fBAM', '-g', genome_size]

        proc = SP.Popen(args, stdout=SP.PIPE, stderr=SP.PIPE)
        proc.wait()
        sout, serr = proc.communicate()
        self.logger.debug(sout)
        self.logger.debug(serr)
        self.logger.info('Finished running MACS for treatment %s.' % (treatment_bam_fn) )
  
if __name__ == '__main__':
    import multiprocessing as MP

    logger = openLogger(log_fn='pipeline.log', log_name='pipeline')

    (options, extra_args) = parseOptions()

    for fn in options.top_data_dir:
        os.path.isdir(fn)
        data_dirs.append(fn)

    def addSample(sampleList, data_dir, cpu_count):
        sample = SeqSample(sample_dir=data_dirs[0], cpu_count, single_str='R1', paired_str='R2')
        sampleList.append(sample)

    sampleList = []
    sample_num = len(data_dirs)
    pool = MP.Pool(num_cpu=sample_num)
    cpu_count = (MP.cpu_count / sample_num)
    pool.map(addSample, [sampleList, data_dir, cpu_count])




    # gunzip and concatenate the fq.gz files
    # fqs is a dictionary with the concatenated filename for 
    # single and paired files
    sample.gunzipFastqs()
    
    sample.trimFastqFiles(keep_untrimmed_fastq=False)

    sample.mapBWA()
    sys.exit()
    
    sample.runMACS(treat_fn=sample.bam_fn, control_fn='', macs_bin='/g/software/bin/macs14')

    
    for bamfn in bamfns:
        makeNormalizedBedGraph(bamfn, species, insertSize=350)
        cleanupSeqFiles(fqs, bamfn, keepsam=False)

