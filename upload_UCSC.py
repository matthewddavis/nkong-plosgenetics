#!/usr/bin/python
import sys
import os
import shutil

def add_header_to_bed(fn, dest, optargs):
    '''.bed files from MACS need the UCSC header line added in order to track
    properly.  This add a minimal header using the file name.
    '''
    print 'Adding header to the bed/bedGraph file %s...' % (fn)

    if not 'samplename' in optargs:
        optargs['samplename'] = fn.split('/')[-1].split('.bed')[0]
    if not 'sampledesc' in optargs:
        optargs['sampledesc'] = ''
    if not 'color' in optargs:
        optargs['color'] = '255,0,0'
    if not 'visibility' in optargs:
        optargs['visibility'] = 'full'

    header = 'track type=%(ftype)s name="%(samplename)s" description="%(sampledesc)s" visibility="%(visibility)s" color=%(color)s autoScale=on maxHeightPixels=100:24:21 graphType=bar alwaysZero=on yLineOnOff=on gridDefault=on windowingFunction=mean\n' % optargs
    print '\t' + header
    infh = open(fn, 'r')
    lines = infh.read()
    if lines[0].startswith('track'):
        print "This file already has a header...not modifying it."
        return fn
    infh.close()
    outfn = dest + optargs['samplename'] + '.' + optargs['ftype']
    print 'Writing updated file to %s' % (outfn)
    outfh = open(outfn, 'w')
    outfh.write(header)
    outfh.write(lines)
    print '...done.'

    return outfn


def upload_file(fn, hgsid, optargs, dest='/home/matt/public_html/tmp/', localserveraddr = 'http://gonzo.qb3.berkeley.edu/~matt/tmp/'):
    '''Copies the given filename to a tmpfile in the http directory, uploads it 
    with wget, and deletes the tmpfile when done.
    '''
    # get the local filename
    if fn.find('/') != -1:
        localfn = fn.split('/')[-1]
    else:
        localfn = fn

    if optargs['ftype'] == 'bed':
        newfn = add_header_to_bed(fn, dest, optargs)
        localurl = localserveraddr + newfn.split('/')[-1]
    elif optargs['ftype'] == 'bedGraph':
        newfn = add_header_to_bed(fn, dest, optargs)
        localurl = localserveraddr + newfn.split('/')[-1]
    else:
        print 'Uploading file as is (no header added)...'
        print 'Copying %s to %s...' % (fn, dest + localfn)
        shutil.copyfile(fn, dest + localfn)
        print '...done.'
        localurl = localserveraddr + localfn
    
    remoteserveraddr = 'http://genome.ucsc.edu/cgi-bin/hgCustom?hgsid=' 
    remoteurl = remoteserveraddr + hgsid + '&hgct_customText=' + localurl

    cmd = "wget --tries=1 --quiet -O /dev/null '" + remoteurl + "'"
    print '#######', remoteurl, cmd
    print "Pushing to \n\t%s \nfrom \n\t%s with command \n\t%s" % (remoteurl, localurl, cmd)
    os.system(cmd)
    print "\nCompleted push."

    #os.remove(dest + localfn)


if __name__ == '__main__':
    
    if (len(sys.argv) < 3):
        print "Usage:  python upload_UCSC.py [filename] [hgsid]; \n\
optional: trackname, samplename, sampledesc, color \
can be assigned with a csv list (i.e. samplename,'Sample A')."
        sys.exit()
    fn = sys.argv[1]
    hgsid = sys.argv[2]
    optargs = dict([arg.split(',', 1) for arg in sys.argv[3:]])
    optargs['ftype'] = fn.split('.')[-1]
    upload_file(fn, hgsid, optargs)


