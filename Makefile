quality_trim:
	for f in *.gz; do \
		gunzip ${f}; \
	done
	qsub trim_6147JAAXX_2_1_pf.sh
	qsub trim_6147JAAXX_3_1_pf.sh
	qsub trim_6147JAAXX_6_1_pf.sh
	qsub trim_6147JAAXX_7_1_pf.sh

velveth:
	qsub velveth_6147JAAXX_2_1.sh
	qsub velveth_6147JAAXX_3_1.sh
	qsub velveth_6147JAAXX_6_1.sh
	qsub velveth_6147JAAXX_7_1.sh

velvetg:
	qsub velvetg_line6i.sh
	qsub velvetg_line6u.sh
	qsub velvetg_line7i.sh
	qsub velvetg_line7u.sh

oases:
	qsub oases_line6i.sh
	qsub oases_line6u.sh
	qsub oases_line7i.sh
	qsub oases_line7u.sh

bowtie_build:
	bowtie2-build -f chick.fa chick3_bowtie2

tophat:
	qsub tophat_6147JAAXX_2_1_pf_trim.sh
	qsub tophat_6147JAAXX_3_1_pf_trim.sh
	qsub tophat_6147JAAXX_6_1_pf_trim.sh
	qsub tophat_6147JAAXX_7_1_pf_trim.sh

install:
	apt-get install -y ncbi-blast+
	wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/blat/blat
	mv blat /usr/local/bin
	git clone https://github.com/likit/gimme.git
	cd gimme; git checkout v.0.97; python setup.py install
	apt-get install -y python-biopython
	apt-get install -y samtools
	wget -O seqclean.tgz http://sourceforge.net/projects/seqclean/files/seqclean-x86_64.tgz/download
	tar xvfz seqclean.tgz
	cd seqclean-x86_64; export PATH=$PATH:$PWD
	cd source; tar xfvz cd-hit-v4.5.4-2011-03-07.tgz; cd cd-hit-v4.5.4-2011-03-07; make && make install
	git clone https://github.com/ctb/screed.git
	cd screed; python setup.py install
	git clone https://github.com/ged-lab/khmer.git
	cd khmer; make && make all
	cd khmer/python; python setup.py build
	cd khmer/python; python setup.py install
	cd source; cp Gimme-paper.ipynb /usr/local/notebooks

clean:
	rm -r gimme
	cd source; rm seqclean.tgz
	cd source; rm -r seqclean-x86_64
	cd source; rm -r cd-hit-v4.5.4-2011-03-07

clean-transcripts:
	cd /mnt/data/data; seqclean se_6u_local.fa
	cd /mnt/data/data; seqclean se_6i_local.fa
	cd /mnt/data/data; seqclean se_7u_local.fa
	cd /mnt/data/data; seqclean se_7i_local.fa
	cd /mnt/data/data; seqclean se_6u_global.fa
	cd /mnt/data/data; seqclean se_6i_global.fa
	cd /mnt/data/data; seqclean se_7u_global.fa
	cd /mnt/data/data; seqclean se_7i_global.fa

remove-redundant-transcripts:
	cd /mnt/data/data; cd-hit-est -d 0 -c 1.0 -M 8000 -i se_6u_local.fa.clean -o se_6u_local.clean.nr.fa
	cd /mnt/data/data; cd-hit-est -d 0 -c 1.0 -M 8000 -i se_6i_local.fa.clean -o se_6i_local.clean.nr.fa
	cd /mnt/data/data; cd-hit-est -d 0 -c 1.0 -M 8000 -i se_7u_local.fa.clean -o se_7u_local.clean.nr.fa
	cd /mnt/data/data; cd-hit-est -d 0 -c 1.0 -M 8000 -i se_7i_local.fa.clean -o se_7i_local.clean.nr.fa
	cd /mnt/data/data; cd-hit-est -d 0 -c 1.0 -M 8000 -i se_6u_global.fa.clean -o se_6u_global.clean.nr.fa
	cd /mnt/data/data; cd-hit-est -d 0 -c 1.0 -M 8000 -i se_6i_global.fa.clean -o se_6i_global.clean.nr.fa
	cd /mnt/data/data; cd-hit-est -d 0 -c 1.0 -M 8000 -i se_7u_global.fa.clean -o se_7u_global.clean.nr.fa
	cd /mnt/data/data; cd-hit-est -d 0 -c 1.0 -M 8000 -i se_7i_global.fa.clean -o se_7i_global.clean.nr.fa

find-unique-transcripts:
	#cd /mnt/data/data; python /mnt/data/gimme/src/utils/assembly-diff-2.py se_6u_local.clean.nr.fa se_6u_global.clean.nr.fa
	cd /mnt/data/data; python /mnt/data/gimme/src/utils/assembly-diff-2.py se_6i_local.clean.nr.fa se_6i_global.clean.nr.fa
	cd /mnt/data/data; python /mnt/data/gimme/src/utils/assembly-diff-2.py se_7u_local.clean.nr.fa se_7u_global.clean.nr.fa
	cd /mnt/data/data; python /mnt/data/gimme/src/utils/assembly-diff-2.py se_7i_local.clean.nr.fa se_7i_global.clean.nr.fa
	cd /mnt/data/data; python /mnt/data/gimme/src/utils/assembly-diff-2.py se_6u_global.clean.nr.fa se_6u_local.clean.nr.fa
	cd /mnt/data/data; python /mnt/data/gimme/src/utils/assembly-diff-2.py se_6i_global.clean.nr.fa se_6i_local.clean.nr.fa
	cd /mnt/data/data; python /mnt/data/gimme/src/utils/assembly-diff-2.py se_7u_global.clean.nr.fa se_7u_local.clean.nr.fa
	cd /mnt/data/data; python /mnt/data/gimme/src/utils/assembly-diff-2.py se_7i_global.clean.nr.fa se_7i_local.clean.nr.fa

preprocess: clean-transcripts remove-redundant-transcripts find-unique-transcripts

global_assembly: quality_trim velveth velvetg oases

local_assembly: bowtie_build tophat

run-blastx:

construct-gene-models:
	cd /mnt/data/data; cat *clean.nr.fa > all_clean.fa
	cd /mnt/data/data; cd-hit-est -T 0 -d 0 -c 1.0 -M 8000 -i all_clean.fa -o all_clean.nr.fa
