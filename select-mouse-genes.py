#!/usr/bin/env python
'''Description'''
import sys

def parse(chrom_list_file):
    chroms = set()
    for line in open(chrom_list_file):
        chroms.add(line.strip())
    return chroms


def main():
    '''Main function'''

    chrom_list_file = sys.argv[1]
    model_file = sys.argv[2]
    chromosomes = parse(chrom_list_file)
    for line in open(model_file):
        chr = line.split('\t')[0]
        if chr in chromosomes:
            print line.strip()


if __name__=='__main__':
    main()

