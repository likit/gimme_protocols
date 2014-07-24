#!/usr/bin/env python
'''Description'''
import sys
from Bio import SeqIO

def count(infile):
    '''Returns a number sequences and a total

    size of all sequences in the input FASTA file.
    '''

    size =  0
    num = 0
    for rec in SeqIO.parse(infile, 'fasta'):
        size += len(rec)
        num += 1

    return num, size


def main():
    '''Main function'''

    infile = sys.argv[1]  # FASTA file
    print 'Total sequences = %d\nTotal size = %d' % \
            count(infile)


if __name__=='__main__':
    main()

