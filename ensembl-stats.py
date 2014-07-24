#!/usr/bin/env python
'''Description'''
import sys


def run(infile):
    '''return basic stats for the gene models'''

    genes = set()
    tranx = set()
    for line in open(infile):
        attrs = line.split('\t')[-1].split(';')
        gene = attrs[0].strip().split()[1].replace('"','')
        trn = attrs[1].strip().split()[1].replace('"','')
        genes.add(gene)
        tranx.add(trn)

    print 'Total genes = ', len(genes)
    print 'Total transcripts = ', len(tranx)

def main():
    '''Main function'''

    infile = sys.argv[1]  # GTF file
    run(infile)


if __name__=='__main__':
    main()

