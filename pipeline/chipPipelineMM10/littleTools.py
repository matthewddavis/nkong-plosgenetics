#!/g/software/bin/python-2.7
import sys

def openNice(fn, attr='r'):
    '''
    Just opens a file nicely.
    '''
    try:
        fh = open(fn, attr)
        return fh
    except:
        print "Error opening file %s." % (fn)

