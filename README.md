##Setup working environment

+ Run

    git clone https://github.com/likit/gimme_protocols.git protocol

###All required packages

+ Blat *
+ BLAST+ 2.2.25 *
+ Bowtie 1.0.0
+ Bowtie2*
+ Tophat2*
+ Cufflinks2*
+ Seqclean *
+ ESTscan 3.0.3*
+ Gimme package 0.97 (+ Networkx, Pysam, Pygr)
+ Python 2.7
+ Biopython
+ Velvet 1.2.03*
+ Oases 0.2.06*
+ CDHIT 4.5.4*
+ Condetri 2.1*
+ Samtools 0.1.18
+ khmer
+ screed
+ RSEM 1.2.7

(*) _These packages are not required to run code in the Ipython notebook._

##Data processing

_Note: Preprocessing may take a long time.
Some steps need to be run on a computer cluster.
You can skip this step and use pre-preprocessed data
to reproduce results in this notebook._

###Quality trimming
Run condetri:

    make -f ~/gimme-protocol/makefile protocol=~/gimme-protocol quality-trim

The quality score cutoff is 30 (default) and the
first 10 nucleotides are hard trimmed.

###Global assembly

Run Velveth and velvetg with kmer=21,23,25,27,29,31:

    make -f ~/gimme-protocol/makefile protocol=~/gimme-protocol velveth-global

do not proceed until the all velveth jobs are finished.

    make -f ~/gimme-protocol/makefile protocol=~/gimme-protocol velvetg-global

Then, clean and remove redundant transcripts by running

    make preprocess

Preprocessing includes the following steps:

+ **clean-transcripts :** clean transcripts with seqclean
+ **remove-redundant-transcripts :** remove redundant transcripts with CD-HIT
+ **find-unique-transcripts :** find unique transcripts between global and local assembly

###EST and mRNAs alignments

Run BLAT to align chicken mRNAs and ESTs:

    make -f ~/gimme-protocol/makefile protcol=~/gimme-protocol run-blat-mrna-est

Sort and select the best alignments (highest score) from mRNAs and ESTs.

    make -f ~/gimme-protocol/makefile protcol=~/gimme-protocol sort-mrna-est-alignments
