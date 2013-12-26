##Setup working environment

To setup an EC2 machine, do the following

###On Amazon EC2
+ Launch EC2 instance with **beacon-2012.09.03 (ami-c17ec8a8)**. Choose **m1.large** for instance type.
+ Create volume with (SNAPSHOT) and mount the volume to /mnt/data/
+ Login as a root
+ Run **git clone https://github.com/likit/gimme_protocols.git protocols**
+ Run **make -f protocols/Makefile install** in **/root** directory
+ Open Ipython notebook by entering https://youramazonpublicdns in your browser (make sure HTTPS rule is included in your security group)
+ Login to the notebook with a password "beacon"

###All required packages

+ Blat *
+ BLAST+ 2.2.25 *
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

(*) _These packages are not required to run code in the Ipython notebook._
