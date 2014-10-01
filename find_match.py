#!/usr/bin/env python
'''Description'''
import sys
from Bio.Blast import NCBIXML

def find_match(xmlfile):
    aligned_bases = 0
    records = NCBIXML.parse(open(xmlfile))
    op = open(xmlfile+'.out', 'w')
    for record in records:
        if record.alignments: # if matched,
            alm = record.alignments[0]
            hsp = alm.hsps[0]
            aligned_bases += hsp.align_length
            if hsp.expect < 1e-20:
                print >> op, '%s\t%s\t%d' % \
                        (record.query,
                            alm.title,
                            hsp.align_length)
    op.close()

def main():
    '''Main function'''

    xmlfile = sys.argv[1]
    find_match(xmlfile)


if __name__=='__main__':
    main()

