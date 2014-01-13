##Setup working environment

To setup an EC2 machine, do the following

###On Amazon EC2
+ Launch EC2 instance with **beacon-2012.09.03 (ami-c17ec8a8)**. Choose **m1.large** for instance type.
+ Create volume with (SNAPSHOT) and mount the volume to /mnt/data/
+ Login as a root and go to **/root**
+ Run **git clone https://github.com/likit/gimme_protocols.git protocols**
+ Run **make -f protocols/Makefile install** in **/root** directory
+ Open Ipython notebook by entering https://youramazonpublicdns in your browser (make sure HTTPS rule is included in your security group)
+ Login to the notebook with a password "beacon"

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

###Data preprocessing

_Note: Preprocessing may take a long time. You can skip this step and use pre-preprocessed data to reproduce results in this notebook._

Run **make preprocess** to preprocess transcripts and other required data to be used in this notebook.

This step will take a while, so we recommend running it in Shell.

In addition, each step can be run separately with its own make command. For example, running **make clean-transcript** will only run Seqclean to clean transcripts.

Preprocessing includes the following steps:

+ **clean-transcripts :** clean transcripts with seqclean
+ **remove-redundant-transcripts :** remove redundant transcripts with CD-HIT
+ **find-unique-transcripts :** find unique transcripts between global and local assembly

These steps are not included in preprocessing. We recommend using pre-preprocessed data:

+ **run-blastx :** run BLASTX against mouse proteins
+ **construct-gene-models :** run Gimme to construct gene models
