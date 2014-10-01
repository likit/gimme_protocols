#!/usr/bin/env python
'''Summarizes the number of canonical and non-canonical splice sites

from splice sites file.
'''
import sys
import re
from pygr import seqdb


pattern = r'(\w+):(\d+)-(\d+)'
pattern = re.compile(pattern)

def get_sp_seq(infile, genome):
    can1 = 0  # GT-AG
    can2 = 0  # GC-AG
    noncano = 0  # non-canonical
    total = 0
    for line in infile:
        ss, strand = line.split()
        chrom, start, end = pattern.findall(ss)[0]
        start = int(start) - 1
        end = int(end) 
        total += 1

        donor = str(genome[chrom][start:start+2])
        acceptor = str(genome[chrom][end-1:end+1])

        if donor == 'GT' and acceptor == 'AG':
            can1 += 1
        elif donor == 'GC' and acceptor == 'AG':
            can2 += 1
        elif donor == 'CT' and acceptor == 'AC':  # reverse strand
            can1 += 1
        elif donor == 'CT' and acceptor == 'TC':  # reverse strand
            can2 += 1
        else:
            noncano += 1
            print ss  # print out non-canonical splice sites

    return can1, can2, noncano, total


def main():
    '''Main function'''

    genome = seqdb.SequenceFileDB(sys.argv[1])  # genome sequence
    try:
        infile = open(sys.argv[2])  # splice sites file from gimme/compare_junction.py
    except IndexError, IOError:
        # no input file given or cannot open the file,
        # read data from stdin instead
        infile = sys.stdin

    can1, can2, noncano, total = get_sp_seq(infile, genome)
    print >> sys.stderr, 'Total splice sites = %d' % total
    print >> sys.stderr, 'Total major canonical sites (GT,AG) = %d (%.2f)' % \
            (can1, float(can1)/total*100)
    print >> sys.stderr, 'Total minor canonical sites (GC,AG) = %d (%.2f)' % \
            (can2, float(can2)/total*100)
    print >> sys.stderr, 'Total non canonical sites = %d (%.2f)' % \
            (noncano, float(noncano)/total*100)



if __name__=='__main__':
    main()

