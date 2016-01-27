#!/usr/bin/python
import sys
import os
import shutil
import logging

def openLogger():
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    # console log is INFO and above
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    # ch.setLevel(logging.INFO)
    # create formatter and add it to the handlers
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    ch.setFormatter(formatter)
    # add the handlers to the logger
    logger.addHandler(ch)

    return logger

def upload_file(dest='/home/matt/public_html/ucsc/', local_server_addr='http://dounce.icmb.utexas.edu/~matt/ucsc/', \
                logger=openLogger(), **kwargs):
    '''
    Copies the bed_fn to a temporary file in the http directory,
    and then uploads it with wget, and deletes the tmpfile when done.

    Accepts either a .bed or .bedGraph file as bed_fn.

    kwargs are:
        bed_fn
        track_name
        uscs_id
        sampl_ename
        sample_desc
        color
        visibility
    '''

    def _getArgs():
        args = {}
        args['bed_fn'] = kwargs.get('bed_fn')
        args['ucsc_id'] = kwargs.get('ucsc_id')
        args['sample_name'] = kwargs.get('sample_name')
        if args['sample_name'] == '':
            args['sample_name'] = args['bed_fn'].split('/')[-1].split('.bed')[0]
        args['sample_desc'] = kwargs.get('sample_desc')
        args['color'] = kwargs.get('color')
        args['visibility'] = kwargs.get('visibility')
        for i in args:
            print args[i]
        return args

    def _add_header_to_bed():
        '''
        bed and bedGraph files need the UCSC header line added in order to track
        properly. This add a minimal header using the file name.
        '''
        logger.info('Adding header to the bed/bedGraph file %s...' % (args['bed_fn']) )

        header = 'track type=%(file_type)s name="%(sample_name)s" description="%(sample_desc)s" visibility="%(visibility)s" \
        color=%(color)s autoScale=on maxHeightPixels=100:24:21 graphType=bar alwaysZero=on yLineOnOff=on gridDefault=on windowingFunction=mean\n' % (args)
        in_fh = open(args['bed_fn'], 'r')
        lines = in_fh.read()
        if lines[0].startswith('track'):
            logger.warning("This %s file already has a header...not modifying it." % (args['bed_fn']) )
        else:
            logger.info('Adding header to %s file.' % (args['file_type']) )
            logger.info('%s' % (header) )
            out_fn = dest + args['sample_name'] + '.' + args['file_type']
            logger.info('Writing updated file to %s' % (out_fn) )
            out_fh = open(out_fn, 'w')
            out_fh.write(header)
            out_fh.write(lines)
            logger.info('...done writing updated file.')

        return out_fn
    args = _getArgs()
    
    if '/' in args['bed_fn']:
        rel_fn = args['bed_fn'].split('/')[-1]
    else:
        rel_fn = args['bed_fn']

    if args['bed_fn'].endswith('.bed'):
        args['file_type'] = 'bed'
        new_fn = _add_header_to_bed()
        local_url = local_server_addr + '/' + new_fn.split('/')[-1]
    elif args['bed_fn'].endswith('.bedGraph'):
        args['file_type'] = 'bedGraph'
        new_fn = _add_header_to_bed()
        local_url = local_server_addr + '/' + new_fn.split('/')[-1]
    else:
        self.logging.info('Uploading file as-is (no header added)...')
        self.logging.info('Copying %s to %s...' % (args['bed_fn'], dest + '/' + local_fn))
        shutil.copyfile(fn, dest + '/' + local_fn)
        self.logging.info('...done.')
        local_url = local_server_addr + local_fn
    
    remote_server_addr = 'http://genome.ucsc.edu/cgi-bin/hgCustom?hgsid='
    remote_url = remote_server_addr + args['ucsc_id'] + '&hgct_customText=' + local_url

    cmd = "wget --tries=1 --quiet -O /dev/null '" + remote_url + "'"
    logger.info('Remote URL: ' + remote_url + cmd)
    logger.info("Pushing to \n\t%s \nfrom \n\t%s with command \n\t%s" % (remote_url, local_url, cmd) )
    os.system(cmd)
    logger.info("Completed push.")

def parseOptions():
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("--ucsc_id", dest="ucsc_id", 
                      help="the UCSC Genome Browser session ID")
    parser.add_option("--bedfile", dest="bed_fn",
                      help="the bed or bedGraph file to loaded into UCSC", metavar="FILE")
    parser.add_option("--trackname", dest="track_name", default='',
                      help="the UCSC track name")
    parser.add_option("--samplename", dest="sample_name", default='',
                      help="the UCSC sample name")
    parser.add_option("--sampledesc", dest="sample_desc", default='',
                      help="the UCSC sample description")
    parser.add_option("--color", dest="color", default='255,0,0',
                      help="the UCSC display color for the track")
    parser.add_option("--visibility", dest="visibility", default='full',
                      help="the UCSC display visibility for the track")
    

    (options, args) = parser.parse_args()
    return (options, args)

if __name__ == '__main__':
    logger = openLogger()
    (options, extra_args) = parseOptions()
    upload_file(**options.__dict__)


