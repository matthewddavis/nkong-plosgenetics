#!/g/software/bin/python-2.7
import sys
import os
from littleTools import openNice

class ChipExoSample():
    '''
    '''
    def __init__(self):
        pass





def gunzipFastqs(sample_dir, read_type='single', species='mm10', out_prefix = 'FASTQ', single_substr='R1', paired_substr='R2'):
    '''
    Unzips and concatenates the FASTQ files from an Illumina experiment
    for downstream analysis.

    Returns a dict of the names of the single and paired end file names for other functions to use.

    The concatenated file is stored in a subdirectory <outdir> of the 
    sample directory <sampledir>.

    The substr identifies the read file as single or paired end.
    '''
    from glob import glob
    import gzip 
    out_dir = sample_dir + '/' + '_'.join([out_prefix, read_type, species])

    print "Uncompressing files for %s to %s..." % (sample_dir, out_dir)
    out_fns = {}
    if os.path.isdir(out_dir):
        sout_fn = glob(out_dir + '/' + '*' + single_substr + '.all.fastq')
        if sout_fn:
            out_fns['single'] = sout_fn[0]
        else:
            return False
        if pout_fn:
            out_fns['paired'] = pout_fn[0]
        return out_fns
        print "#####"
        print "Warning: It looks like the unzipped FASTQ files are already present in %s" % (out_dir)
        print "#####"
        pout_fn = glob(out_dir + '/' + '*' + paired_substr + '.all.fastq')


    else:
        os.mkdir(out_dir)
    
    fastq_gzs = glob(sample_dir + '/' + '*.gz')

    samplename = fastq_gzs[0].split('/')[-1].split(single_substr)[0] 
    
    single_fq_gzs, paired_fq_gzs = [], []
    for fastq_gz in fastq_gzs:
        if single_substr in fastq_gz:
            single_fq_gzs.append(fastq_gz)
        elif paired_substr in fastq_gz:
            paired_fq_gzs.append(paired_substr)
    for fastq_gz in fastq_gzs:
        if (fastq_gz not in single_fq_gzs) and (fastq_gz not in paired_fq_gzs):
            print "#####"
            print "WARNING:  Skipping %s because it doesn't contain a single or paired end FASTQ substring!"            
            print "#####"
    
    out_fns = {}

    if len(single_fq_gzs) > 0:
        sout_fn = out_dir + '/' + sample_name + single_substr + '.all.fastq'
        os.system('zcat ' + ' '.join(single_fq_gzs) + '> ' + sout_fn)
        out_fns['single'] = sout_fn
    if len(paired_fq_gzs) > 0:
        pout_fn = out_dir + '/' + sample_name + paired_substr + '.all.fastq'
        os.system('zcat ' + ' '.join(paired_fq_gzs) + '> ' + pout_fn)
        out_fns['paired'] = pout_fn

    return out_fns
