PACKAGES = install_blast install_blat install_gimme install_biopython \
		install_samtools install_seqclean install_cdhit install_khmer_screed \
		install_tophat2 install_bowtie2 install_velvet install_oases \
		install_condetri

.PHONY: $(PACKAGES)

preinstall:

	mkdir /mnt/source

install: clean preinstall $(PACKAGES)

	cd /usr/local/notebooks; ln -sf /root/protocol/notebooks.ipynb

install_blast:

	apt-get install -y ncbi-blast+

install_blat:

	cd /mnt/source; wget http://genome-test.cse.ucsc.edu/~kent/exe/linux/blatSuite.zip; \
		unzip blatSuite.zip

install_gimme:

	cd /mnt/source; git clone https://github.com/likit/gimme.git
	cd /mnt/source/gimme; git checkout v.0.97; python setup.py install

install_biopython:

	apt-get install -y python-biopython

install_samtools:

	apt-get install -y samtools

install_seqclean:

	cd /mnt/source; wget -O seqclean.tgz http://sourceforge.net/projects/seqclean/files/seqclean-x86_64.tgz/download
	cd /mnt/source; tar xvfz seqclean.tgz;

install_cdhit:

	cd /mnt/source; wget https://cdhit.googlecode.com/files/cd-hit-v4.5.4-2011-03-07.tgz
	cd /mnt/source/cd-hit-v4.5.4-2011-03-07; make install

install_khmer_screed:

	cd /mnt/source; tar xfvz cd-hit-v4.5.4-2011-03-07.tgz; cd cd-hit-v4.5.4-2011-03-07; make && make install
	cd /mnt/source; git clone https://github.com/ctb/screed.git; cd screed; python setup.py install
	cd /mnt/source; git clone https://github.com/ged-lab/khmer.git; cd khmer; make
	cd /mnt/source/khmer/python; python setup.py install

install_velvet:

	cd /mnt/source; wget http://www.ebi.ac.uk/~zerbino/velvet/velvet_1.2.03.tgz; tar xvfz velvet_1.2.03.tgz
	cd /mnt/source/velvet_1.2.03/; make 'MAXKMERLENGTH=57'
	cd /mnt/source/velvet_1.2.03/; cp shuffleSequences_fast*.pl velveth velvetg /usr/local/bin

install_oases:

	cd /mnt/source; wget http://www.ebi.ac.uk/~zerbino/oases/oases_0.2.06.tgz; tar xvfz oases_0.2.06.tgz
	cd /mnt/source; cd oases_0.2.06; make 'VELVET_DIR=/mnt/source/velvet_1.2.03' 'MAXKMERLENGTH=57'
	cd /mnt/source/oases_0.2.06; cp oases /usr/local/bin/

install_tophat2:

	cd /mnt/source; wget http://tophat.cbcb.umd.edu/downloads/tophat-2.0.5.Linux_x86_64.tar.gz
	cd /mnt/source; tar xvfz tophat-2.0.5.Linux_x86_64.tar.gz
	cd /mnt/source/tophat-2.0.5.Linux_x86_64; find * -executable -exec cp '{}' /usr/local/bin \;
	
install_bowtie2:

	cd /mnt/source; wget http://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.0.0-beta5/bowtie2-2.0.0-beta5-linux-x86_64.zip
	cd /mnt/source; unzip bowtie2-2.0.0-beta5-linux-x86_64.zip
	cd /mnt/source/bowtie2-2.0.0-beta5; cp bowtie2 bowtie2-align bowtie2-build bowtie2-inspect /usr/local/bin
	
install_condetri:

	cd /mnt/source; wget https://condetri.googlecode.com/files/condetri_v2.1.pl

install_cufflinks:

	cd /mnt/source; wget http://cufflinks.cbcb.umd.edu/downloads/cufflinks-2.0.0.Linux_x86_64.tar.gz
	cd /mnt/source; tar xvfz cufflinks-2.0.0.Linux_x86_64.tar.gz
	cd /mnt/source; mv cufflinks-2.0.0.Linux_x86_64 cufflinks2

clean:

	rm -r /mnt/source
