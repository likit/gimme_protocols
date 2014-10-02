##Setup working environment

Run:

    git clone https://github.com/likit/gimme_protocols.git protocol
    export PROTOCOL=<path to protocol>
    export GIMMEDIR=<path to gimme>
    alias runmake="make -f $PROTOCOL/makefile protocol=$PROTOCOL gimmedir=$GIMMEDIR"

###All required packages

+ Blat
+ BLAST+ 2.2.25
+ Bowtie 1.0.0
+ Tophat 1.3.1
+ Cufflinks 2.0.0
+ Seqclean
+ ESTscan 3.0.3
+ Gimme package 0.97 (need a new release)
+ Python 2.7
+ Biopython
+ Pysam
+ Velvet 1.2.03
+ Oases 0.2.06
+ CDHIT 4.5.6
+ Condetri 2.1
+ Samtools 0.1.18
+ khmer 0.8
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

    runmake quality-trim

The quality score cutoff is 30 (default) and the
first 10 nucleotides are hard trimmed.

###Global assembly

Run Velveth and Velvetg with kmer=21,23,25,27,29,31:

    runmake velveth-global

Do not proceed until the velveth job is finished.

    runmake velvetg-global

Run Oases:

    runmake oases-global

Then, clean and remove redundant transcripts by running

Run Velveth and Velvetg for OasesM:

    runmake velveth-M

Do not proceed until the velveth job is finished.

    runmake velvetg-M

Run OasesM:

    runmake oases-M

Clean oases-M transcripts and align them to the chicken genome:

    runmake run-blat-oases-M

Sort and filter alignments:

    runmake sort-oases-M-alignments

###Local assembly

Map reads to chicken genome:

    runmake tophat

Extract reads from chromosomes:

    runmake extract-reads

Run Velveth local:

    runmake local-velveth

Run Velvetg local:

    runmake local-velvetg

Run Oases local:

    runmake local-oases

Combine and clean transcripts:

    runmake combine-global-assembly-transcripts combine-local-assembly-transcripts
    runmake clean-transcripts

Remove redundant transcripts:

    runmake remove-redundant-all-assembly

Align all transcripts to the chicken genome

    runmake run-blat-all-assembly

Remove redundant global and local transcripts:

    runmake remove-redundant-global-assembly remove-redundant-local-assembly

Align global and local transcripts to the chicken genome

    runmake run-blat-global-assembly run-blat-local-assembly run-blat-all-assembly

###Build assembly gene models

Construct gene models from global assembly:

    runmake construct-gene-models-global

Construct gene models from local assembly:

    runmake construct-gene-models-local

Construct gene models from global and local assembly:

    runmake construct-gene-models-global-local

###Find unique transcripts from global and local assembly

    runmake find-unique-transcripts
    runmake get-long-unique-regions

###Build Cufflinks models

Run cufflinks:

    runmake cufflinks

Run cufflinks merge:

    runmake cufflinks-merge

Construct global+local+cufflinks models:

    runmake construct-gene-models-global-local-cufflinks


###Mouse data

Quality trim:

Interleave reads:

    runmake interleave-mouse-reads

Run Velveth and Velvetg:

    runmake velveth-mouse

Note, mouse data is strand-specific.

Do not proceed until the velveth job is finished.

    runmake velvetg-mouse

Run Oases:

    runmake oases-mouse

Run Tophat:

    runmake tophat-mouse

Extract reads:

    runmake extract-reads-mouse

Run mouse local assembly:

    runmake local-velveth-mouse
    runmake local-velvetg-mouse
    runmake local-oases-mouse

Combine, clean and remove redundant transcripts:

    runmake combine-global-mouse-assembly-transcripts
    runmake combine-local-mouse-assembly-transcripts
    runmake clean-remove-redundant-mouse-transcripts

Align mouse assembly to the mouse genome:

    runmake run-blat-all-mouse-assembly

Note, mouse dataset is strand-specific.
This version of Gimme does not support strand-specific data internally.
Therefore, models of positive- and negative-strand genes are built separately.

Construct mouse Gimme models:

    runmake construct-gene-models-mouse

###EST and mRNAs alignments

Run BLAT to align chicken mRNAs and ESTs:

    runmake run-blat-mrna-est

Sort and select the best alignments (highest score) from mRNAs and ESTs.

    runmake sort-mrna-est-alignments

###Reads mapping statistics

Build bowtie index for gene models.

    runmake build-bowtie-index-gimme-models
    runmake build-bowtie-index-cufflinks-models
    runmake gimme-models-map
    runmake cufflinks-map

###Splice junctions analysis

    runmake run-blat-global-two-kmers
    runmake sort-blat-two-kmers
    runmake find-splice-junctions

    runmake count-spliced-reads-ensembl
    runmake count-spliced-reads-gimme
    runmake count-spliced-reads-cufflinks

###Homology analysis

    runmake create-mouse-db
    runmake gimme-vs-mouse-blastx
    runmake gimme-vs-mouse-blastp
    runmake run-blastx-unique-regions
    runmake find-mouse-match

    runmake translate-cufflinks
    runmake blast-cufflinks-vs-mouse
    runmake parse-cufflinks-vs-mouse-blast

    runmake get-sequences-cufflinks-gimme-models
    runmake translate-cufflinks-gimme
    runmake blast-cufflinks-gimme-vs-mouse
    runmake parse-cufflinks-gimme-vs-mouse-blast

###Error profiles

    runmake reads-error-profile
