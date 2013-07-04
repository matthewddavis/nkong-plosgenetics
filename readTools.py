#!/g/software/bin/python-2.7
import sys

def gunzipFastqs(sample_dir, read_type, species='mm10', out_prefix = 'FASTQ', single_substr='R1', paired_substr='R2'):
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
