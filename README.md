##Setup working environment

+ Run

    git clone https://github.com/likit/gimme_protocols.git protocol

###All required packages

+ Blat
+ BLAST+ 2.2.25
+ Bowtie 1.0.0
+ Tophat 1.3.1
+ Cufflinks2 2.0.0
+ Seqclean
+ ESTscan 3.0.3
+ Gimme package 0.97 (need a new release)
+ Python 2.7
+ Biopython
+ Velvet 1.2.03
+ Oases 0.2.06
+ CDHIT 4.5.4
+ Condetri 2.1
+ Samtools 0.1.18
+ khmer
+ screed
+ RSEM 1.2.7

(*) _These packages are not required to run code in the Ipython notebook._

##Raw data processing

_Note: These steps take days to finish and
some of them need to be run on a big-mem computer.
You can skip this step and use pre-preprocessed data
to reproduce results and plots included in the paper._

###Quality trimming
Run condetri:

    make -f $PROTOCOL/makefile protocol=$PROTOCOL quality-trim

The quality score cutoff is 30 (default) and the
first 10 nucleotides are hard trimmed.

###Global assembly

Run Velveth and Velvetg with kmer=21,23,25,27,29,31:

    make -f $PROTOCOL/makefile protocol=$PROTOCOL velveth-global

Do not proceed until the velveth job is finished.

    make -f $PROTOCOL/makefile protocol=$PROTOCOL velvetg-global

Run Oases:

    make -f $PROTOCOL/makefile protocol=$PROTOCOL oases-global

Then, clean and remove redundant transcripts by running

    make preprocess

Preprocessing includes the following steps:

+ **clean-transcripts :** clean transcripts with seqclean
+ **remove-redundant-transcripts :** remove redundant transcripts with CD-HIT
+ **find-unique-transcripts :** find unique transcripts between global and local assembly

Run Velveth and Velvetg for OasesM:

    make -f $PROTOCOL/makefile protocol=$PROTOCOL velveth-M

Do not proceed until the velveth job is finished.

    make -f $PROTOCOL/makefile protocol=$PROTOCOL velvetg-M

Run OasesM:

    make -f $PROTOCOL/makefile protocol=$PROTOCOL oases-M

Clean oases-M transcripts and align them to the chicken genome:

    make -f $PROTOCOL/makefile protocol=$PROTOCOL run-blat-oases-M

###Local assembly

Map reads to chicken genome:

    make -f $PROTOCOL/makefile protocol=$PROTOCOL tophat

Extract reads from chromosomes:

    make -f $PROTOCOL/makefile protocol=$PROTOCOL extract-reads

Run Velveth local:

    make -f $PROTOCOL/makefile protocol=$PROTOCOL local-velveth

Run Velvetg local:

    make -f $PROTOCOL/makefile protocol=$PROTOCOL local-velvetg

Run Oases local:

    make -f $PROTOCOL/makefile protocol=$PROTOCOL local-oases

Combine and clean transcripts:

    make -f $PROTOCOL/makefile protocol=$PROTOCOL gimmedir=$GIMMEDIR  combine-global-assembly-transcripts
    make -f $PROTOCOL/makefile protocol=$PROTOCOL gimmedir=$GIMMEDIR  combine-local-assembly-transcripts
    make -f $PROTOCOL/makefile clean-transcripts

Remove redundant transcripts:

    make -f $PROTOCOL/makefile protocol=$PROTOCOL remove-redundant-all-assembly

Align all transcripts to the chicken genome

    make -f $PROTOCOL/makefile protocol=$PROTOCOL run-blat-all-assembly

Remove redundant global and local transcripts:

    make -f $PROTOCOL/makefile protocol=$PROTOCOL remove-redundant-global-assembly
    make -f $PROTOCOL/makefile protocol=$PROTOCOL remove-redundant-local-assembly

Align global and local transcripts to the chicken genome

    make -f $PROTOCOL/makefile protocol=$PROTOCOL run-blat-global-assembly
    make -f $PROTOCOL/makefile protocol=$PROTOCOL run-blat-local-assembly

Construct gene models from global assembly:

    make -f $PROTOCOL/makefile protocol=$PROTOCOL gimmedir=$GIMMEDIR construct-gene-models-global

Construct gene models from local assembly:

    make -f $PROTOCOL/makefile protocol=$PROTOCOL gimmedir=$GIMMEDIR construct-gene-models-local

Construct gene models from global and local assembly:

    make -f $PROTOCOL/makefile protocol=$PROTOCOL gimmedir=$GIMMEDIR construct-gene-models-global-local

###Build Cufflinks models

Run cufflinks:

    make -f $PROTOCOL/makefile protocol=$PROTOCOL cufflinks

Run cufflinks merge:

    make -f $PROTOCOL/makefile cufflinks-merge

###Mouse data

Quality trim:

Interleave reads:

    make -f $PROTOCOL/makefile protocol=$PROTOCOL interleave-mouse-reads

Run Velveth and Velvetg:

    make -f $PROTOCOL/makefile protocol=$PROTOCOL velveth-mouse

Note, mouse data is strand-specific.

Do not proceed until the velveth job is finished.

    make -f $PROTOCOL/makefile protocol=$PROTOCOL velvetg-mouse

Run Oases:

    make -f $PROTOCOL/makefile protocol=$PROTOCOL oases-mouse

Run Tophat:

    make -f $PROTOCOL/makefile protocol=$PROTOCOL tophat-mouse

Extract reads:

    make -f $PROTOCOL/makefile protocol=$PROTOCOL extract-reads-mouse

Run mouse local assembly:

    make -f $PROTOCOL/makefile local-velveth-mouse
    make -f $PROTOCOL/makefile local-velvetg-mouse
    make -f $PROTOCOL/makefile local-oases-mouse

Combine, clean and remove redundant transcripts:

    make -f $PROTOCOL/makefile gimmedir=$GIMMEDIR combine-global-mouse-assembly-transcripts
    make -f $PROTOCOL/makefile gimmedir=$GIMMEDIR combine-local-mouse-assembly-transcripts
    make -f $PROTOCOL/makefile protocol=$PROTOCOL clean-remove-redundant-mouse-transcripts

Align mouse assembly to the mouse genome:

    make -f $PROTOCOL/makefile protocol=$PROTOCOL run-blat-all-mouse-assembly

###EST and mRNAs alignments

Run BLAT to align chicken mRNAs and ESTs:

    make -f $PROTOCOL/makefile protcol=$PROTOCOL run-blat-mrna-est

Sort and select the best alignments (highest score) from mRNAs and ESTs.

    make -f $PROTOCOL/makefile protcol=$PROTOCOL sort-mrna-est-alignments

###Splice junctions analysis

Align global transcripts k21 and k31 to chicken genome:

    make -f $PROTOCOL/makefile protcol=$PROTOCOL run-blat-global-two-kmers

###Reads mapping statistics

Build bowtie index for gene models.

