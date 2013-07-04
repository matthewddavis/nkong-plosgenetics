#!/usr/bin/python
import sys


def trimFastq(fastqfn, adapterseq):
    '''
    '''

    fastqfh = open(fastqfn,'r')

    outfh = open(fastqfn.split('.')[0] + '.adaptercleaned.fastq', 'w')

    count = 0
    while (True):
        for i in range(4):
            seqname = fastqfh.readline().strip()
            if seqname == '': return
            read = removeAdapter(fastqfh.readline().strip(), adapterseq)
            if read == '':
                read = 'N'
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

def main():
    if len(sys.argv) < 2:
        print "Usage:  trimReads.py <infiles>"
        sys.exit()
    print sys.argv[1:]
    trimFastq(sys.argv[1], sys.argv[2])


if __name__ == '__main__':
    main()

