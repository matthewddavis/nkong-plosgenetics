#!/g/software/bin/python-2.7
import sys
import logging
import os
import inspect
import subprocess
import threading


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
    # console log is INFO and above
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
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

def cleanupSeqFiles(sample, keep_unzipped_fastq=False, keep_sam=False, keep_bigwig=False, keep_bed=False,\
                    keep_wig=False, keep_tmp_bedGraph=False, keep_bedGraph=True):
    '''
    Deletes intermediate files of large size.
    '''
    fns = []
    if not keep_unzipped_fastq:
        unzipped_fq_fns = glob(sample.sample_dir + '/*.fastq')
        fns.extend(unzipped_fq_fns)

    if not keep_sam:
        sai_fns = glob(sample.bwa_dir + '/*.sai')
        fns.extend(sai_fns)
        sam_fns = glob(sample.bwa_dir + '/*.sam')
        fns.extend(sam_fns)
        
    if not keep_bigwig:
        bw_fns = glob(sample.bwa_dir + '/*.bw')
        fns.extend(bw_fns)

    if not keep_bed:
        bed_fns = glob(sample.bwa_dir + '/*.bed')
        fns.extend(bed_fns)
        
    if not keep_wig:
        wig_fns = glob(sample.bwa_dir + '/*.wig')
        fns.extend(wig_fns)

    if not keep_tmp_bedGraph:
        tmp_bedGraph_fns = glob(sample.macs_dir + '/*.bedGraph.tmp')
        fns.extend(tmp_bedGraph_fns)

    if not keep_bedGraph:
        bedGraph_fns = glob(sample.bwa_dir + '/*.bedGraph')
        fns.extend(bedGraph_fns)

    for fn in fns:
        sample.logger.info("Removing file %s" % (fn))
        try:
            os.remove(fn)
        except OSError as e:
            sample.logger.error("ERROR: couldn't remove file %s due to %s" % (fn, e) )

class SeqSample(object):
    '''
    '''
    def __init__(self, sample_dir, cpu_count=1, single_str='R1', paired_str='R2', remove_duplicate_reads=True):
        self.sample_dir = sample_dir
        self.cpu_count = str(cpu_count)
        self.single_str = single_str
        self.paired_str = paired_str
        self.remove_duplicate_reads = remove_duplicate_reads

        self.attrs = {}
        self._loadSampleAttrs()

        self.use_untrimmed = False

        try:
            self.logger = openLogger(log_fn=sample_dir + '/pipeline.log', log_name=self.attrs['sample_id'])
        except KeyError:
            pdb.set_trace()

        self._checkSampleAttrs()
        self._setControlDataFile()

        self.bwa_idx = ''
        self.bwa_dir = ''

    def __getstate__(self):
        d = dict(self.__dict__)
        del d['logger']
        return d

    def _setControlDataFile(self):
        if self.attrs['control_sample'] == '':
            self.control_data_bam_fn = ''
        else:
            control_dir = 'data/' + self.attrs['control_sample']
            if self.remove_duplicate_reads:
                control_bam_fn = glob(control_dir + '/bwa/*rmdup.bam')[0]
            else:
                control_bam_fn = glob(control_dir + '/bwa/*.bam')[0]
                self.control_data_bam_fn = control_bam_fn
            assert os.path.isfile(control_bam_fn)
            self.control_data_bam_fn = control_bam_fn
            
    def _loadSampleAttrs(self):
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
        insert_size; 150

        This function parses these attrs into the SeqSample object.
        '''
        sample_info_fn = self.sample_dir + '/' + 'sample_info'
        for line in open(sample_info_fn):
            if line.startswith('#'):
                continue
            line = line.split(';')
            line[1] = line[1].strip()
            if line[0] == 'sample_id':
                line[1] = line[1].replace(' ', '_')
            self.attrs[line[0]] = line[1]

    def _checkSampleAttrs(self,):
        self.logger.info('Sample loaded:')
        for attr in self.attrs:
            self.logger.info('\t%s : %s' % (attr, self.attrs[attr]))
            
        for required_key in ['sample_id', 'common_name', 'data_type', 'read_type', 'read_length', 'adapter_seq', 'species', 'insert_size']:
            if required_key not in self.attrs:
                self.logger.error('Required info %s is missing from sample_info file.' % (required_key))
                sys.exit()

    def gunzipFastqs(self):
        '''
        Unzips and concatenates the FASTQ files from an Illumina experiment
        for downstream analysis.

        The concatenated file is stored in a subdirectory <outdir> of the
        sample directory <sampledir>.

        The single_str or paired_str identifies the read file as single or paired end.
        For VCGSL reads, this is R1 for single and R2 for paired.
        '''
        def _checkForFastq(out_dir):
            '''
            Check to see if the uncompressed files already exist from a previous run.

            The out_dir can be anywhere, but usually is the sample_dir.
            '''
            read_type = self.attrs['read_type']
            out_fns = glob(out_dir + '/' + '*' + 'all.fastq')
            # out_fns = glob(out_dir + '/' + '*' + '.fastq')
            if len(out_fns) == 0:
                return False

            # There is at least one all.fastq, but there should be at most two all.fastq files
            # (in the case of paired-end reads)
            assert len(out_fns) < 3
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
            self.logger.info("We already have uncompressed files at %s" % (self.sample_dir) )
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
                    #self.logged.debug('Found the adapter at pos %s in read %s' % (idx, read))
                    return read[0:idx]
                else:
                    idx = read.rfind(adapter_seq[0:3])
                    if idx != -1:
                        pad = len(read) - idx
                        idx = read.rfind(adapter_seq[0:pad])
                        if idx != -1:
                            #self.logger.debug('Found the incomplete adapter at pos %s in read %s' % (idx, read) )
                            return read[0:idx]
                        else:
                            #self.logger.debug('Not even incomplete adapter found in read %s' % (read) )
                            return read
                    else:
                        #self.logger.debug('No adapter found in read %s' % (read) )
                        pass
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
            #single_fastq_fh = open(self.single_fastq_fn, 'r')
            single_trimmed_fastq_fn = self.single_fastq_fn.replace('.all.', '.trimmed.')
            if _checkForTrimmedFastq():
                self.single_trimmed_fastq_fn = single_trimmed_fastq_fn
                self.paired_trimmed_fastq_fn = None
                return
            else:
                _trimFastqFile(self.single_fastq_fn)
                self.paired_trimmed_fastq_fn = None
                
        elif self.attrs['read_type'] == 'paired':
            #single_fastq_fh = open(self.single_fastq_fn, 'r')
            single_trimmed_fastq_fn = self.single_fastq_fn.replace('.all.', '.trimmed.')
            #paired_fastq_fh = open(self.paired_fastq_fn, 'r')
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
            sys.exit()

        if not keep_untrimmed_fastq:
            _removeUntrimmedFastqs()


    def checkForBAMFiles(self,):
       '''
       Check for any BAM files, and if we find them, skip generation.
       '''
       bam_fns = glob(self.sample_dir + '/bwa/*.bam')
       if len(bam_fns) > 0:
           self.logger.warn("WARNING: We already have .bam files such as %s, so we're skipping read mapping." % (bam_fns[0]) )
           return True


    def indexReadsBWA(self, bwa_bin='bin/bwa'):
        '''
        Index the reads with BWA.
        '''
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

        def _checkForSAI():
            sai_fns = glob(self.bwa_dir + '/*.sai')
            if sai_fns == []:
                self.logger.info('No .sai files found for %s...creating .sai file(s) now' % (self.attrs['sample_id']) )
                return False
            else:
                self.logger.warn('Warning: Found some previous .sai files found for %s' % (self.attrs['sample_id']) )
                if self.attrs['read_type'] == 'single':
                    for sai_fn in sai_fns:
                        if self.single_str in sai_fn:
                            self.logger.warn('Using %s for the single-end .sai file.' % (sai_fn) )
                            self.single_sai_fn = sai_fn
                    if self.single_sai_fn != None:
                        return True
                elif self.attrs['read_type'] == 'paired':
                    for sai_fn in sai_fns:
                        if self.single_str in sai_fn:
                            self.logger.warn('Using %s for the single-end .sai file.' % (sai_fn) )
                            self.single_sai_fn = sai_fn
                        if self.paired_str in sai_fn:
                            self.logger.warn('Using %s for the paired-end .sai file.' % (sai_fn) )
                            self.paired_sai_fn = sai_fn
                    if (self.single_sai_fn != None) and (self.paired_sai_fn != None):
                        return True

        _createBWADir()
            
        bwa_dir = self.bwa_dir

        if self.use_untrimmed:
            fq_fns = [self.single_fastq_fn, self.paired_fastq_fn]
            self.logger.info("Using untrimmed FASTQ files for BWA indexing.")
        else:
            fq_fns = [self.single_trimmed_fastq_fn, self.paired_trimmed_fastq_fn]
            self.logger.info("Using trimmed FASTQ files for BWA indexing.")
            
        for fq_fn in fq_fns:
            if fq_fn == None:
                continue
            else:
                self.logger.info("Making BWA index file for fastq file %s (sample %s)..." % (fq_fn, self.attrs['sample_id']) )
            # make a directory to store things in
            sai_fn = bwa_dir + '/' + fq_fn.replace(self.sample_dir + '/', '').replace('.fastq', '.sai')
            # if the file is already there, skip to the next file
            #if os.path.isfile(sai_fn):
            #    continue
            if self.single_str in sai_fn:
                self.single_sai_fn = sai_fn
            elif self.paired_str in sai_fn:
                self.paired_sai_fn = sai_fn
            else:
                self.logger.error('!!! The read_type is not present or not interpretable in the sample_info file.')
                self.logger.error('!!! Exiting on line %s...' % (inspect.currentframe().f_lineno) )
            
            if _checkForSAI():
                return

            sai_log_fn = sai_fn.replace('.fastq', '.sai.log')
            sai_log_fh = open(sai_log_fn, 'w')
            args = [bwa_bin, 'aln', '-t', str(self.cpu_count), self.bwa_idx, fq_fn, '-f', sai_fn]
            # proc = subprocess.Popen(args, stdout=sai_log_fh, stderr=sai_log_fh)
            proc = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout, stderr = proc.communicate()
            self.logger.info(stdout)
            self.logger.info(stderr)
            sai_log_fh.close()
            self.logger.info("Finished making BWA idx file %s for fastq file %s (%s)." % (sai_fn, fq_fn, self.attrs['sample_id']) )

    def alignBWA(self, bwa_bin='bin/bwa'):
        '''
        Align reads with BWA.
        '''
        def _checkForSam():
            self.logger.info("Checking for a previous .sam file for %s." % (self.attrs['sample_id']) )
            if os.path.isfile(sam_fn):
                self.logger.info("Found a previous .sam file for %s." % (self.attrs['sample_id']) )
                self.sam_fn = sam_fn
                return True

        bwa_idx = self.bwa_idx
        #bwa_dir = self.bwa_dir

        if self.use_untrimmed:
            fqs = [self.single_fastq_fn, self.paired_fastq_fn]
        else:
            fqs = [self.single_trimmed_fastq_fn, self.paired_trimmed_fastq_fn]

        if self.attrs['read_type'] == 'single':
            self.logger.info('Aligning single-end reads for sample %s' % (self.attrs['sample_id']) )
            sam_fn = self.single_sai_fn.replace('.sai', '.sam')
            if _checkForSam():
                return
            #if self.checkForBAMFiles():
            #    return
            sam_log_fn = self.single_sai_fn.replace('.sai', '.sam.log')
            sam_log_fh = open(sam_log_fn, 'w')
            self.logger.info("Aligning single-end reads for sample %s with BWA to %s." % (self.attrs['sample_id'], self.single_sai_fn) )
            args = [bwa_bin, 'samse', bwa_idx, self.single_sai_fn, fqs[0], '-f', sam_fn]
            proc = subprocess.Popen(args, stdout=sam_log_fh, stderr=sam_log_fh)
            #proc.wait()
            proc.communicate()

        elif self.attrs['read_type'] == 'paired':
            # make sure we have a paired fastq file
            assert fqs[1] != None
            self.logger.info('Aligning paired-end reads for sample %s' % (self.attrs['sample_id']) )
            # doesn't matter which sai file you pull the name from,
            # as long as the both exist (and this implicitly checks for them)
            sam_fn = self.paired_sai_fn.replace('.sai', '.sam')
            if _checkForSam():
                return
            #if self.checkForBAMFiles():
            #    return
            sam_log_fn = self.paired_sai_fn.replace('.sai', '.sam.log')
            sam_log_fh = open(sam_log_fn, 'w')
            self.logger.info("Aligning paired-end reads for sample %s with BWA to %s." % (self.attrs['sample_id'], self.paired_sai_fn) )
            args = [bwa_bin, 'sampe', bwa_idx, self.single_sai_fn, self.paired_sai_fn, fqs[0], fqs[1], '-f', sam_fn]
            proc = subprocess.Popen(args, stdout=sam_log_fh, stderr=sam_log_fh)
            proc.communicate()

        self.sam_fn = sam_fn

    def parseSAMByStrand(self):
        '''
        For data types such as ChIP-Exo, we may want to split the reads into two sam files,
        one for each strand. This method does this using the bit flags in the sam files.

        This is typically just done for visualization in a genome browser, but could also
        be used for a custom analysis of strand biased enrichment.
        '''

        def _checkForStrandSamFiles():
            if os.path.isfile(split_sam_watson_fn) and os.path.isfile(split_sam_crick_fn):
                return True

        split_sam_watson_fn = self.sam_fn.replace('.sam', '.watson.sam')
        split_sam_crick_fn = self.sam_fn.replace('sam', 'crick.sam')

        if _checkForStrandSamFiles():
            self.split_sam_watson_fn = split_sam_watson_fn
            self.split_sam_crick_fn = split_sam_crick_fn
            return

        split_sam_watson_fh = open(split_sam_watson_fn, 'w')
        split_sam_crick_fh = open(split_sam_crick_fn, 'w')
        sam_fh = open(self.sam_fn, 'r')

        self.logger.info("Parsing the sam file for sample %s into separate strands." % (self.attrs['sample_id']) )
        line_counter = 0
        for line in sam_fh:
            line_counter += 1
            if (line_counter != 0) and (line_counter % 100000 == 0):
                self.logger.info("So far we've parsed %s lines by strand." % (line_counter) )
            # the sam header lines start with @
            if line.startswith('@'):
                split_sam_watson_fh.write(line)
                split_sam_crick_fh.write(line)
                continue
            # if the read is on the reverse strand, demarcated by the bitwise 16 operation,
            # send it to crick, otherwise, send it to watson
            line = line.split('\t')
            if (int(line[1]) & 16):
                split_sam_crick_fh.write('\t'.join(line))
            else:
                split_sam_watson_fh.write('\t'.join(line))

        self.split_sam_watson_fn = split_sam_watson_fn
        self.split_sam_crick_fn = split_sam_crick_fn

    def parseSAMBySpecies(self):
        '''
        For mixed-species chromatin samples, we split them into two separate .sam files,
        one for each species.
        '''
        # other are reads that map with more than one best alignment
        #    and we toss these, because we want to avoid any possible
        #    cross-mapping

        def _checkForParsedSAMFiles():
            if os.path.isfile(self.dm3_sam_fn):
                self.logger.warn("We already have species-parsed .sam files, such as %s." % (self.dm3_sam_fn) )
                return True
            else:
                self.logger.info("Parsing sam file into three species (dm3, dp3, other)")
                self.logger.info("New sam files are: %s\n%s\n%s\n" % (self.dm3_sam_fn, self.dp3_sam_fn, self.other_sam_fn) )
                return False
            
        self.dm3_sam_fn = self.sam_fn.replace('.sam', '.dm3.sam')
        self.dp3_sam_fn = self.sam_fn.replace('.sam', '.dp3.sam')
        self.other_sam_fn = self.sam_fn.replace('.sam', '.other.sam')

        if _checkForParsedSAMFiles():
            return

        dm3_sam_fh = open(self.dm3_sam_fn, 'w')
        dp3_sam_fh = open(self.dp3_sam_fn, 'w')
        other_sam_fh = open(self.other_sam_fn, 'w')
        sam_fh = open(self.sam_fn, 'r')

        counter = 0
        dm3_counter = 0
        dp3_counter = 0
        other_counter = 0
        for line in sam_fh:
            if (counter % 1e5 == 0) and (counter > 0):
                self.logger.info("Parsed %s reads by species so far..." % (counter) )
            # first deal with the sam header lines
            if line.startswith('@'):
                if 'DMEL' in line:
                    dm3_sam_fh.write(line.replace('DMEL', ''))
                    continue
                elif 'DPSE' in line:
                    dp3_sam_fh.write(line.replace('DPSE', ''))
                    continue
            # get the lines for each species' chromosomes
            # and take only the ones with only one best alignment (X0:i:1)
            counter += 1
            if ('\tDMEL' in line) and ('X0:i:1\t' in line):
                dm3_counter += 1
                dm3_sam_fh.write(line.replace('\tDMEL', '\t'))
            elif ('\tDPSE' in line) and ('X0:i:1\t' in line):
                dp3_counter += 1
                dp3_sam_fh.write(line.replace('\tDPSE', '\t'))
            else:
                other_counter += 1
                other_sam_fh.write(line)

        self.logger.info("%s reads sent to %s" % (dm3_counter, self.dm3_sam_fn) )
        self.logger.info("%s reads sent to %s" % (dp3_counter, self.dp3_sam_fn) )
        self.logger.info("%s reads sent to %s" % (other_counter, self.other_sam_fn) )
        self.logger.info("%s total reads parsed by species" % (counter) )

    def samToBam(self, sam_fn):
        '''
        Converts sam to bam and returns the bam_fn so that any sam_fn can
        be converted and stored externally.
        '''
        def _sortBAM(bam_fn):
            '''
            Sorts a bam file.
            '''
            self.logger.info('Sorting bam file %s...' % (bam_fn) )
            # the replace() below is required b/c samtools sort wants to always add .bam to the file
            args = ['samtools', 'sort', bam_fn, bam_fn.replace('.bam', '')]
            proc = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            proc.wait()
            sout, serr = proc.communicate()
            self.logger.debug(sout)
            self.logger.debug(serr)
            self.logger.info('Finished sorting bam file %s.' % (bam_fn) )

        def _checkForBam(bam_fn):
            if os.path.isfile(bam_fn):
                self.logger.warn("We already have a bam file this sam file. \n\tbam: %s\n\tsam: %s" % (bam_fn, sam_fn) )
                return True
            else:
                self.logger.info('Converting %s to .bam format...' % (sam_fn) )
                return False
            
        ## SamToBam main body begins here
        bam_fn = sam_fn.replace('.sam', '.bam')
        if (_checkForBam(bam_fn)):
            return bam_fn
        args = ['samtools', 'view', '-bS', sam_fn, '-o', bam_fn]
        proc = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        proc.wait()
        sout, serr = proc.communicate()
        self.logger.debug(sout)
        self.logger.debug(serr)
        self.logger.info('Finished converting %s to .bam format...' % (sam_fn) )

        _sortBAM(bam_fn)
        self.makeBAI(bam_fn)

        return bam_fn

    def makeBAI(self, bam_fn):
        '''
        Creates the bam idx file required to read the .bam file.
        '''
        bai_fn = bam_fn.replace('.bam', '.bai')
        self.logger.info('Making .bai file for %s...' % (self.attrs['sample_id']) )
        args = ['samtools', 'index', bam_fn, bai_fn]
        proc = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        proc.wait()
        sout, serr = proc.communicate()
        self.logger.debug(sout)
        self.logger.debug(serr)
        self.logger.info('Finish making .bai file for %s.' % (bam_fn) )

    def removeDuplicateReads(self, bam_fn):
        '''
        In some cases, such as ChIP-Seq, we may want to remove exact
        duplicate reads from the samples.
        '''
        def _checkForDuplicatesRemovedFiles(duplicates_removed_bam_fn):
            if os.path.isfile(duplicates_removed_bam_fn):
                self.logger.warn("We already have bam files with the duplicate reads removed, such as %s" % (duplicates_removed_bam_fn) )
                return True
            else:
                return False

        duplicates_removed_bam_fn = bam_fn.replace('.bam', '.rmdup.bam')
        if (_checkForDuplicatesRemovedFiles(duplicates_removed_bam_fn)):
            return duplicates_removed_bam_fn
        self.logger.info("Removing exact mapped duplicates for sample %s from bam file %s..." % (self.attrs['sample_id'], bam_fn) )
        args = ['samtools', 'rmdup', '-s', bam_fn, duplicates_removed_bam_fn]
        proc = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        proc.wait()
        sout, serr = proc.communicate()
        self.logger.debug(sout)
        self.logger.debug(serr)
        self.logger.info('Finished removing exact mapped duplicates for sample %s from bam file %s' % (self.attrs['sample_id'], bam_fn) )

        # need to remake the .bai file after removing duplicate reads
        self.logger.info('Re-making the .bai file for %s after removing duplicate reads' % (self.attrs['sample_id']) )
        self.makeBAI(duplicates_removed_bam_fn)

        return duplicates_removed_bam_fn

    def mapBWA(self, use_untrimmed=False):
        '''
        Mapping reads with BWA (as opposed to some other read mapper).
        '''
        def _setUntrimmedFlag():
            '''
            If use_untrimmed is set True, the untrimmed FASTQ files will be mapped.
            '''
            if use_untrimmed:
                assert os.path.isfile(self.single_fastq_fn)
            self.use_untrimmed = use_untrimmed

        def _getBWAIDX_fn():
            if self.attrs['species'] == 'dm3':
                bwa_idx = 'reference/bwaIndices/dm3/dm3.fa'
                self.bwa_idx = bwa_idx
            elif self.attrs['species'] == 'dp3':
                bwa_idx = 'reference/bwaIndices/dp3/dp3.fa'
                self.bwa_idx = bwa_idx
            elif self.attrs['species'] == 'dm3dp3':
                bwa_idx = 'reference/bwaIndices/dm3dp3/dm3dp3Merge.fa'
                self.bwa_idx = bwa_idx
            elif self.attrs['species'] == 'mm10':
                bwa_idx = 'reference/bwaIndices/mm10_repeat_masked/mm10.masked.fa'
                self.bwa_idx = bwa_idx
            else:
                print self.attrs
                self.logger.error("Error: invalid species name for %s" % (self.attrs['species']) )
                self.logger.error('!!! Exiting on line %s...' % (inspect.currentframe().f_lineno) )
                sys.exit()

        _setUntrimmedFlag()
        _getBWAIDX_fn()
        bwa_idx = self.bwa_idx

        # map and align the reads
        self.indexReadsBWA()
        self.alignBWA()

        # make the BAM files
        if self.attrs['data_type'] == 'chip-exo':
            self.parseSAMByStrand()
            self.bam_fn = self.samToBam(self.sam_fn)
            self.split_bam_waston_fn = self.samToBam(self.split_sam_watson_fn)
            self.split_bam_crick_fn = self.samToBam(self.split_sam_crick_fn)
        elif self.attrs['data_type'] == 'chip-seq':
            if self.attrs['species'] == 'dm3dp3':
                self.parseSAMBySpecies()
                self.dm3_bam_fn = self.samToBam(self.dm3_sam_fn)
                self.dp3_bam_fn = self.samToBam(self.dp3_sam_fn)
                self.other_bam_fn = self.samToBam(self.other_sam_fn)
                # set the macs_treatment_bam_fn first to this bam_fn,
                # if we want to use the duplicates_removed file instead,
                # it gets set below
                self.dm3_macs_treatment_bam_fn = self.dm3_bam_fn
                self.dp3_macs_treatment_bam_fn = self.dp3_bam_fn
                self.other_macs_treatment_bam_fn = self.other_bam_fn
                if self.remove_duplicate_reads:
                    self.dm3_duplicates_removed_bam_fn = self.removeDuplicateReads(self.dm3_bam_fn)
                    self.dp3_duplicates_removed_bam_fn = self.removeDuplicateReads(self.dp3_bam_fn)
                    self.other_duplicates_removed_bam_fn = self.removeDuplicateReads(self.other_bam_fn)
                    self.dm3_macs_treatment_bam_fn = self.dm3_duplicates_removed_bam_fn
                    self.dp3_macs_treatment_bam_fn = self.dp3_duplicates_removed_bam_fn
                    self.other_macs_treatment_bam_fn = self.other_duplicates_removed_bam_fn
            elif self.attrs['species'] == 'dm3':
                self.dm3_bam_fn = self.samToBam(self.sam_fn)
                self.dm3_macs_treatment_bam_fn = self.dm3_bam_fn
                if self.remove_duplicate_reads:
                    self.dm3_duplicates_removed_bam_fn = self.removeDuplicateReads(self.dm3_bam_fn)
                    self.dm3_macs_treatment_bam_fn = self.dm3_duplicates_removed_bam_fn

    def runMACS(self, treatment_bam_fn, macs_bin='/usr/bin/macs14'):
        '''
        Uses MACS to call peaks in a dataset.

        The treatment_bam_fn can must be a sorted bam file.

        If control_sample_path is unspecified, MACS14 will use the no input control model
        '''
        #def _getControlBamFn():
        #    '''
        #    Gets the control bam filename for the control sample path specified in the
        #    sample_info file.
        #    '''
        #    bam_fns = glob(control_sample_path + '/' + self.bwa_dir + '/' + '*.bam')
        #    # need to add some code here for the split read bam files (15 July 2013)
        #    assert len(bam_fns) == 1
        #    return bam_fns[0]

        def _checkForMACSFile():
            if os.path.isfile(macs_fn):
                self.logger.warning("Warning: we already have MACS output at %s" % macs_fn)
                return True
            else:
                return False
        

        # make sure MACS is installed
        assert os.path.isfile(macs_bin)

        if self.attrs['species'] == 'mm10':
            genome_size = '1.87e9'
        elif self.attrs['species'] == 'dm3':
            genome_size = '1.2e8'
        elif self.attrs['species'] == 'dm3dp3':
            genome_size = '1.2e8'
        else:
            self.logger.error('!!! The species is not specified or not interpretable in the sample_info file.')
            self.logger.error('!!! Exiting on line %s...' % (inspect.currentframe().f_lineno) )
            sys.exit()

        self.macs_dir = self.bwa_dir.replace('/bwa', '/macs')
        if not os.path.isdir(self.macs_dir):
            os.mkdir(self.macs_dir)

        macs_fn = treatment_bam_fn.replace('.bam', '.macs')
        macs_fn = macs_fn.split('/')[-1]
        macs_fn = self.macs_dir + '/' + macs_fn
        if _checkForMACSFile():
            self.macs_peaks_bed_fn = macs_fn + '_peaks.bed'
            return

        if self.control_data_bam_fn == '':
            self.logger.info("Running MACS with %s as the treatment and no control file..." % (treatment_bam_fn) )
            args = [macs_bin, '-t', treatment_bam_fn, '-n', macs_fn, '-fBAM', '-g', genome_size]
        else:
            # control_bam_fn = _getControlBamFn()
            print "\tRunning MACS14 with %s as the treatment and %s as the control..." % (treatment_bam_fn, self.control_data_bam_fn)
            args = [macs_bin, '-t', treatment_bam_fn, '-c', self.control_data_bam_fn, '-n', macs_fn, '-fBAM', '-g', genome_size]

        macs_log_fh = open(macs_fn, 'w')
        #proc = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        proc = subprocess.Popen(args, stdout=macs_log_fh, stderr=macs_log_fh)
        #proc.wait()
        #import pdb
        #pdb.set_trace()
        sout, serr = proc.communicate()
        #self.logger.debug(sout)
        #self.logger.debug(serr)
        self.logger.info('Finished running MACS for treatment %s.' % (treatment_bam_fn) )
        self.macs_peaks_bed_fn = macs_fn + '_peaks.bed'

def splitMACSPeaks(sample, peak_splitter_bin='bin/PeakSplitter'):
    '''
    '''
    def _formatBedGraph():
        '''
        Need to modify the bedGraph file to have a header
        for PeakSplitter to work.
        '''
        bg_fh = open(bg_fn, 'r')
        first_line = bg_fh.readline()
        if not first_line.startswith('track'):
            new_bg_fn = bg_fn + '.tmp'
            new_bg_fh = open(new_bg_fn, 'w')
            new_bg_fh.write('track type=bedGraph\n')
            while True:
                line = bg_fh.readline()
                if line == '':
                    break
                else:
                    new_bg_fh.write(line)
            bg_fh.close()
            new_bg_fh.close()
            return new_bg_fn

    def _removeSplitMACSHeader():
        fh = open(sample.split_peaks_bed_fn, 'r')
        lines = fh.readlines()
        lines.pop(0)
        fh.close()
        fh = open(sample.split_peaks_bed_fn, 'w')
        fh.writelines(lines)

    peaks_bed_fn = sample.macs_peaks_bed_fn
    bg_fn = sample.bg_fn
    macs_dir = sample.macs_dir
    
    bg_fn = _formatBedGraph()

    bg_fh = open(bg_fn, 'r')
    first_line = bg_fh.readline()

    sample.logger.info("Splitting peaks for %s..." % (peaks_bed_fn) )

    subpeaks_dir = macs_dir + '/peaksplitter'
    if not os.path.isdir(subpeaks_dir):
        os.mkdir(subpeaks_dir)

    args=[peak_splitter_bin, '-p', peaks_bed_fn, '-w', bg_fn, '-o', subpeaks_dir, '-f']
    proc = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    (stdout, stderr) = proc.communicate()
    sample.logger.info(stdout)
    sample.logger.info(stderr)
    sample.split_peaks_bed_fn = subpeaks_dir + '/' + peaks_bed_fn.split('/')[-1].replace('.trimmed', '.subpeaks.trimmed')

    _removeSplitMACSHeader()
    
def makeScaledBedGraph(seq_sample, bam_fn, scaled_lib_size=10e7):
    '''
    '''
    def _checkForBedGraphFile():
        bg_fn = bam_fn.replace('.bam', '.bedGraph')
        if os.path.isfile(bg_fn):
            seq_sample.bg_fn = bg_fn
            return True

    def _loadChromSizes():
        if seq_sample.attrs['species'] == 'mm10':
            chrom_size_fn = 'reference/mm10.chrom.sizes'
            assert os.path.isfile(chrom_size_fn)
            seq_sample.logger.info('Using chromosome sizes from %s' % (chrom_size_fn) )
        elif seq_sample.attrs['species'] == 'dm3':
            chrom_size_fn = 'reference/dm3.chrom.sizes'
            assert os.path.isfile(chrom_size_fn)
            seq_sample.logger.info('Using chromosome sizes from %s' % (chrom_size_fn) )
        elif seq_sample.attrs['species'] == 'dm3dp3':
            chrom_size_fn = 'reference/dm3.chrom.sizes'
            assert os.path.isfile(chrom_size_fn)
            seq_sample.logger.info('Using chromosome sizes from %s' % (chrom_size_fn) )
        else:
            seq_sample.logger.error('ERROR: cannot find a valid chromosome size file.')
            sys.exit()

        chrom_sizes = {}
        chrom_size_fh = open(chrom_size_fn, 'r')
        for line in chrom_size_fh:
            line = line.strip().split('\t')
            if line == ['']:
                continue
            else:
                seq_sample.logger.info(line)
                chrom_sizes[line[0]] = int(line[1])
        return chrom_sizes, chrom_size_fn

    def _calculateMappedLibrarySize():
        '''
        Use samtools to calculate the mapped library size.
        '''
        seq_sample.logger.info('Calculating the mapped library size for %s' % (bam_fn) )
        args = ['samtools', 'idxstats', bam_fn]
        proc = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        # releveant output is from stderr
        stderr, stdout = proc.communicate()
        size = 0
        for line in stderr.split('\n'):
            line = line.split()
            if line:
                size += int(line[2])
        seq_sample.logger.info('A total of %s are mapped in %s' % (size, bam_fn) )

        with open(bam_fn.replace('.bam', '.mappedLibrarySize'), 'w') as size_fh:
            size_fh.write(str(size))

        return size

    def _bamToBed():
        '''
        Use UCSC tools to write bed file from bam file.
        '''
        bed_fn = bam_fn.replace('.bam', '.bed')
        bed_fh = open(bed_fn, 'w')
        args = ['bamToBed', '-i', bam_fn]
        proc = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        bedPIPE = proc.communicate()[0].split('\n')

        count = 0
        insert_size = int(seq_sample.attrs['insert_size'])

        for line in bedPIPE:
            count += 1
            line = line.split()
            # somewhere the bamToBed program introduces a blank line, this handles it
            try:
                chrom = line[0]
            except IndexError:
                continue
            start = int(line[1])
            stop = int(line[2])
            strand = line[5]
            read_len = stop - start
            if strand == '+':
                stop = start + insert_size
                if stop > chrom_sizes[chrom]:
                    stop = chrom_sizes[chrom]
            elif strand == '-':
                stop = stop + read_len
                start = stop - insert_size
                if start < 1:
                    start = 1
                if stop > chrom_sizes[chrom]:
                    stop = chrom_sizes[chrom]
            bed_fh.write('\t'.join([chrom,str(start),str(stop),strand]) + '\n')
        bed_fh.close()
        return bed_fn

    def _sortBedFile():
        '''
        Use GNU sort to sort the bed file.
        '''
        seq_sample.logger.info("Sorting the bed file %s with GNU sort..." % (bed_fn) )
        args = ['sort', '-k1,1', '-k2,2n', bed_fn, '-o', bed_fn]
        proc = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        sout, serr = proc.communicate()
        if serr:
            seq_sample.logger.debug(serr)

    def _makeWigFile():
        '''
        Make the density coverage file (wig) file.
        '''
        wig_fn = bam_fn.replace('.bam', '.wig')
        wig_fh = open(wig_fn, 'w')

        seq_sample.logger.info("Making the density coverage file with genomeCoverageBed...")
        args = ['genomeCoverageBed', '-i', bed_fn, '-g', chrom_size_fn, '-dz']
        proc = subprocess.Popen(args, stdout=wig_fh, stderr=wig_fh)
        proc.wait()
        wig_fh.close()

        return wig_fn

    def _makeScaledWigFile():
        '''
        Writes scaled density file
        '''
        seq_sample.logger.info('Scaling the wig file.')
        wig_fh = open(wig_fn, 'r')
        scaled_wig_fn = wig_fn.replace('.wig', '.scaled.wig')
        scaled_wig_fh = open(scaled_wig_fn, 'w')

        chrom = ''
        #last_line = ''
        for line in wig_fh:
            line = line.split('\t')
            try:
                cur_chrom = line[0]
                # wig files are 1-based
                pos = int(line[1]) + 1
                count = float(line[2]) * scaled_lib_size / library_size
            except Exception:
                continue
            if not cur_chrom == chrom:
                seq_sample.logger.info("\t...scaling %s..." % (cur_chrom) )
                chrom = cur_chrom
                scaled_wig_fh.write('variableStep chrom=%s\n' % (chrom))
            out_line = '\t'.join([str(pos), str(count)])
            scaled_wig_fh.write(out_line + '\n')
            #last_line = line
        scaled_wig_fh.close()
        wig_fh.close()

        return scaled_wig_fn

    def _makeBigWigFile():
        '''
        Create BigWig
        '''
        bw_fn = bam_fn.replace('.bam', '.bw')
        seq_sample.logger.info("Creating BigWig file %s..." % (bw_fn) )
        args = ['wigToBigWig', scaled_wig_fn, chrom_size_fn, bw_fn]
        proc = subprocess.Popen(args)
        proc.wait()

        return bw_fn

    def _makeBedGraphFile():
        '''
        Create BedGraph
        '''
        bg_fn = bam_fn.replace('.bam', '.bedGraph')
        seq_sample.logger.info("Creating BedGraph file %s..." % (bg_fn) )
        args = ['bigWigToBedGraph', bw_fn, bg_fn]
        proc = subprocess.Popen(args)
        proc.wait()

        return bg_fn


    if _checkForBedGraphFile():
        return
    seq_sample.logger.info('Creating BedGraph files scaled to a total of %s reads.' % (scaled_lib_size) )
    chrom_sizes, chrom_size_fn = _loadChromSizes()
    library_size = _calculateMappedLibrarySize()
    bed_fn = _bamToBed()
    _sortBedFile()
    wig_fn = _makeWigFile()
    scaled_wig_fn = _makeScaledWigFile()
    bw_fn = _makeBigWigFile()
    bg_fn = _makeBedGraphFile()
    seq_sample.bg_fn = bg_fn

def sendToUCSC(sample, send_bedGraph=True, send_macs_bed=True, send_split_macs_bed=True):
    import upload_UCSC as UCSC
    def _getUCSCID(ucsc_id_fn='./ucsc_id'):
        fh = open(ucsc_id_fn)
        return fh.readline().strip()

    ucsc_id = _getUCSCID()

    if send_bedGraph:
        bedGraph_args = { 'bed_fn' : sample.bg_fn,
                          'track_name' : sample.attrs['common_name'],
                          'ucsc_id' : ucsc_id,
                          'sample_name' : sample.attrs['sample_id'] + '_Scaled_Density',
                          'sample_desc' : sample.attrs['common_name'] + ' Scaled Density',
                          'color' : sample.attrs['ucsc_color'],
                          'visibility' : sample.attrs['ucsc_visibility'],
                          }
        sample.logger.info('Sending read density .bedGraph to UCSC Broswer')
        UCSC.upload_file(logger=sample.logger, **bedGraph_args)
        sample.logger.info('Done sending read density .bedGraph to UCSC Broswer')

    if send_macs_bed:
        bed_args = { 'bed_fn' : sample.macs_peaks_bed_fn,
                     'track_name' : sample.attrs['common_name'],
                     'ucsc_id' : ucsc_id,
                     'sample_name' : sample.attrs['sample_id'] + '_MACS_peaks',
                     'sample_desc' : sample.attrs['common_name'] + ' MACS peaks',
                     'color' : sample.attrs['ucsc_color'],
                     'visibility' : sample.attrs['ucsc_visibility'],
                     }
        sample.logger.info('Sending MACS peaks .bed to UCSC Broswer')
        UCSC.upload_file(logger=sample.logger, **bed_args)
        sample.logger.info('Done sending MACS peaks .bed to UCSC Broswer')

    if send_split_macs_bed:
        split_bed_args = { 'bed_fn' : sample.split_peaks_bed_fn,
                           'track_name' : sample.attrs['common_name'],
                           'ucsc_id' : ucsc_id,
                           'sample_name' : sample.attrs['sample_id'] + '_Split_MACS_peaks',
                           'sample_desc' : sample.attrs['common_name'] + ' Split MACS peaks',
                           'color' : sample.attrs['ucsc_color'],
                           'visibility' : sample.attrs['ucsc_visibility'],
                           }
        sample.logger.info('Sending Split MACS peaks .bed to UCSC Broswer')
        UCSC.upload_file(logger=sample.logger, **split_bed_args)
        sample.logger.info('Done sending Split MACS peaks .bed to UCSC Broswer')
        
def processTreatmentSample(data_dir):
    args = { 'cpu_count' : 16,
             'single_str' : 'R1',
             'paired_str' : 'R2',
             'remove_duplicate_reads' : False,
             }
    sample = SeqSample(data_dir, **args)
    sample.gunzipFastqs()
    sample.trimFastqFiles(keep_untrimmed_fastq=True)
    sample.mapBWA()
    sample.runMACS(treatment_bam_fn=sample.bam_fn, macs_bin='/g/software/bin/macs14')
    makeScaledBedGraph(sample, sample.bam_fn, scaled_lib_size=10e7)
    splitMACSPeaks(sample)
    sendToUCSC(sample, send_bedGraph=True, send_macs_bed=True, send_split_macs_bed=True)
    cleanupSeqFiles(sample, keep_unzipped_fastq=False, keep_sam=False, keep_bigwig=False, keep_bed=False,\
                    keep_wig=False, keep_tmp_bedGraph=False, keep_bedGraph=True)

def processInputSample(data_dir):
    kwargs = { 'cpu_count' : 20,
               'single_str' : 'R1',
               'paired_str' : 'R2',
               'remove_duplicate_reads' : False,
               }

    sample = SeqSample(data_dir, **kwargs)
    #sample.gunzipFastqs()
    #sample.trimFastqFiles(keep_untrimmed_fastq=True)
    #sample.mapBWA()
    #makeScaledBedGraph(sample, sample.bam_fn, scaled_lib_size=10e7)
    #sendToUCSC(sample, send_bedGraph=True, send_macs_bed=False, send_split_macs_bed=False)
    # these don't have bedGraphs, that's why we're "keeping" them
    cleanupSeqFiles(sample, keep_unzipped_fastq=False, keep_sam=False, keep_bigwig=False, keep_bed=False,\
                    keep_wig=False, keep_tmp_bedGraph=True, keep_bedGraph=True)

class ThreadedProcessInputSample(threading.Thread):
    '''
    '''
    def __init__(self, queue):
        threading.Thread.__init__(self)
        self.queue = queue

    def run(self):
        processInputSample(data_dir=self.queue.get())
        self.queue.task_done()

class ThreadedProcessTreatmentSample(threading.Thread):
    '''
    '''
    def __init__(self, queue):
        threading.Thread.__init__(self)
        self.queue = queue

    def run(self):
        processTreatmentSample(data_dir=self.queue.get())
        self.queue.task_done()
        
if __name__ == '__main__':
    import Queue
    import pdb
    from glob import glob
    
    (options, extra_args) = parseOptions()

    treatment_dirs = []
    input_dirs = []
    for fn in glob(options.top_data_dir + '/*'):
        if os.path.isdir(fn):
            if 'Input' in fn:
                input_dirs.append(fn)
            elif 'ChIP_September' in fn:
                treatment_dirs.append(fn)
                
    # processInputSample(input_dirs[1])

    
    queue = Queue.Queue()
    for i in range(2):
        # t = ThreadedProcessInputSample(queue)
        t = ThreadedProcessTreatmentSample(queue)
        t.setDaemon(True)
        t.start()
    # for data_dir in input_dirs:
    for data_dir in treatment_dirs:
        queue.put(data_dir)
    queue.join()
    
    




