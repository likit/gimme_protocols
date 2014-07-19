quality_trim:

	for f in *.gz; do \
		gunzip $$f; \
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

tophat:

	qsub tophat_6147JAAXX_2_1_pf_trim.sh
	qsub tophat_6147JAAXX_3_1_pf_trim.sh
	qsub tophat_6147JAAXX_6_1_pf_trim.sh
	qsub tophat_6147JAAXX_7_1_pf_trim.sh

########################
# On EC2 Machine
########################

tophat_ec2:

	cd /mnt/data; tophat -p 2 -o line6u_tophat chick3_bowtie2 6147JAAXX_2_1_pf_trim.fastq
	cd /mnt/data; tophat -p 2 -o line6i_tophat chick3_bowtie2 6147JAAXX_3_1_pf_trim.fastq
	cd /mnt/data; tophat -p 2 -o line7u_tophat chick3_bowtie2 6147JAAXX_6_1_pf_trim.fastq
	cd /mnt/data; tophat -p 2 -o line7i_tophat chick3_bowtie2 6147JAAXX_7_1_pf_trim.fastq

index_samfiles:

	cd /mnt/data/line6u_tophat; samtools index accepted_hits.bam
	cd /mnt/data/line6i_tophat; samtools index accepted_hits.bam
	cd /mnt/data/line7u_tophat; samtools index accepted_hits.bam
	cd /mnt/data/line7i_tophat; samtools index accepted_hits.bam

extract_reads:

	cd /mnt/data; \
	for dir in line*tophat; do \
		cd $$dir; \
		printf "working on %s:\n" $$dir; \
		for chr in `cat /mnt/chromosomes.list`; do \
			printf "\textracting %s...\n" $$chr; \
			samtools view -b -o $$chr.bam accepted_hits.bam $$chr; \
		done; \
		cd /mnt; \
	done

local_velveth:

	cd /mnt/data; \
	for dir in line*tophat; do \
		cd $$dir; \
		for chr in chr*bam; do \
			velveth `basename $$chr .bam`_asm 21,33,2 -short -bam $$chr; \
		done; \
		cd /mnt; \
	done

local_velvetg:

	cd /mnt/data; \
		for dir in line*tophat; do \
		cd $$dir; \
		for chr in chr*asm*; do \
			velvetg $$chr -read_trkg yes -unused_reads yes; \
		done; \
		cd /mnt; \
	done

local_oases:

	cd /mnt/data; \
		for dir in line*tophat; do \
		cd $$dir; \
		for chr in chr*asm*; do \
			oases $$chr -unused_reads yes; \
		done; \
		cd /mnt; \
	done

combine-transcripts:

	cd /mnt/data/line6u_tophat; \
	for d in chr*asm*; \
		do python /mnt/source/gimme/src/utils/rename_fasta.py $$d/transcripts.fa line6u_local_$$d >> ../line6u_local.fa; \
	done
	cd /mnt/data/line6i_tophat; \
	for d in chr*asm*; \
		do python /mnt/source/gimme/src/utils/rename_fasta.py $$d/transcripts.fa line6i_local_$$d >> ../line6i_local.fa; \
	done
	cd /mnt/data/line7u_tophat; \
	for d in chr*asm*; \
		do python /mnt/source/gimme/src/utils/rename_fasta.py $$d/transcripts.fa line7u_local_$$d >> ../line7u_local.fa; \
	done
	cd /mnt/data/line7i_tophat; \
	for d in chr*asm*; \
		do python /mnt/source/gimme/src/utils/rename_fasta.py $$d/transcripts.fa line7i_local_$$d >> ../line7i_local.fa; \
	done

	cd /mnt/data/; \
		for f in line6u_global_*.transcripts.fa; \
		do python /mnt/source/gimme/src/utils/rename_fasta.py $$f $$(basename $$f .transcripts.fa) >> line6u_global.fa; \
	done
	cd /mnt/data/; \
		for f in line6i_global_*.transcripts.fa; \
		do python /mnt/source/gimme/src/utils/rename_fasta.py $$f $$(basename $$f .transcripts.fa) >> line6i_global.fa; \
	done
	cd /mnt/data/; \
		for f in line7u_global_*.transcripts.fa; \
		do python /mnt/source/gimme/src/utils/rename_fasta.py $$f $$(basename $$f .transcripts.fa) >> line7u_global.fa; \
	done
	cd /mnt/data/; \
		for f in line7i_global_*.transcripts.fa; \
		do python /mnt/source/gimme/src/utils/rename_fasta.py $$f $$(basename $$f .transcripts.fa) >> line7i_global.fa; \
	done

clean-transcripts:

	cd /mnt/data/; /mnt/source/seqclean-x86_64/seqclean line6u_local.fa
	cd /mnt/data/; /mnt/source/seqclean-x86_64/seqclean line6i_local.fa
	cd /mnt/data/; /mnt/source/seqclean-x86_64/seqclean line7u_local.fa
	cd /mnt/data/; /mnt/source/seqclean-x86_64/seqclean line7i_local.fa
	cd /mnt/data/; /mnt/source/seqclean-x86_64/seqclean line6u_global.fa
	cd /mnt/data/; /mnt/source/seqclean-x86_64/seqclean line6i_global.fa
	cd /mnt/data/; /mnt/source/seqclean-x86_64/seqclean line7u_global.fa
	cd /mnt/data/; /mnt/source/seqclean-x86_64/seqclean line7i_global.fa

remove-redundant-transcripts:

	cd /mnt/data; /mnt/source/cd-hit-est -d 0 -c 1.0 -M 8000 -i line6u_local.fa.clean -o line6u_local.fa.clean.nr
	cd /mnt/data; /mnt/source/cd-hit-est -d 0 -c 1.0 -M 8000 -i line6i_local.fa.clean -o line6i_local.fa.clean.nr
	cd /mnt/data; /mnt/source/cd-hit-est -d 0 -c 1.0 -M 8000 -i line7u_local.fa.clean -o line7u_local.fa.clean.nr
	cd /mnt/data; /mnt/source/cd-hit-est -d 0 -c 1.0 -M 8000 -i line7i_local.fa.clean -o line7i_local.fa.clean.nr

	cd /mnt/data; /mnt/source/cd-hit-est -d 0 -c 1.0 -M 8000 -i line6u_global.fa.clean -o line6u_global.fa.clean.nr
	cd /mnt/data; /mnt/source/cd-hit-est -d 0 -c 1.0 -M 8000 -i line6i_global.fa.clean -o line6i_global.fa.clean.nr
	cd /mnt/data; /mnt/source/cd-hit-est -d 0 -c 1.0 -M 8000 -i line7u_global.fa.clean -o line7u_global.fa.clean.nr
	cd /mnt/data; /mnt/source/cd-hit-est -d 0 -c 1.0 -M 8000 -i line7i_global.fa.clean -o line7i_global.fa.clean.nr

preprocess: clean-transcripts remove-redundant-transcripts find-unique-transcripts

global_assembly: quality_trim velveth velvetg oases

local_assembly: tophat index_samfiles

construct-gene-models-global:

	cd /mnt/data; cat *global*clean.nr > all.global.fa.clean
	cd /mnt/data; cd-hit-est -T 0 -d 0 -c 1.0 -M 8000 -i all.global.fa.clean -o all.global.fa.clean.nr
	cd /mnt/data; blat -noHead -out=psl -mask=lower -extendThroughN chick_3.2bit all.global.fa.clean.nr all.global.fa.clean.nr.psl
	cd /mnt/data; sort -k 10 all.global.fa.clean.nr.psl > all.global.fa.clean.nr.psl.sorted
	cd /mnt/data; ../source/pslReps -nohead -singleHit all.global.fa.clean.nr.psl.sorted all.global.fa.clean.nr.psl.best info
	cd /mnt/data; python ../source/gimme/src/gimme.py all.global.fa.clean.nr.psl.best > all.global.fa.clean.nr.bed

construct-gene-models-local:

	cd /mnt/data; cat *local*clean.nr > all.local.fa.clean
	cd /mnt/data; cd-hit-est -T 0 -d 0 -c 1.0 -M 8000 -i all.local.fa.clean -o all.local.fa.clean.nr
	cd /mnt/data; blat -noHead -out=psl -mask=lower -extendThroughN chick_3.2bit all.local.fa.clean.nr all.local.fa.clean.nr.psl
	cd /mnt/data; sort -k 10 all.local.fa.clean.nr.psl > all.local.fa.clean.nr.psl.sorted
	cd /mnt/data; ../source/pslReps -nohead -singleHit all.local.fa.clean.nr.psl.sorted all.local.fa.clean.nr.psl.best info
	cd /mnt/data; python ../source/gimme/src/gimme.py all.local.fa.clean.nr.psl.best > all.local.fa.clean.nr.bed

construct-gene-models-global-local:

	cd /mnt/data; cat *clean.nr > all.fa.clean
	cd /mnt/data; cd-hit-est -T 0 -d 0 -c 1.0 -M 8000 -i all.fa.clean -o all.fa.clean.nr
	cd /mnt/data; blat -noHead -out=psl -mask=lower -extendThroughN chick_3.2bit all.fa.clean.nr all.fa.clean.nr.psl
	cd /mnt/data; sort -k 10 all.fa.clean.nr.psl > all.fa.clean.nr.psl.sorted
	cd /mnt/data; ../source/pslReps -nohead -singleHit all.fa.clean.nr.psl.sorted all.fa.clean.nr.psl.best info
	cd /mnt/data; python ../source/gimme/src/gimme.py all.fa.clean.nr.psl.best > all.fa.clean.nr.bed

clean-up-gene-models-global:

	cd /mnt/data; python /mnt/source/gimme/src/utils/get_transcript_seq.py all.global.fa.clean.nr.bed chick.fa > all.global.fa.clean.nr.bed.fa
	cd /mnt/data; cd-hit-est -T 0 -d 0 -c 0.99 -M 8000 -i all.global.fa.clean.nr.bed.fa -o all.global.fa.clean.nr.bed.fa.nr99
	python ../source/gimme/src/utils/cdhit_transcript.py all.global.fa.clean.nr.bed all.global.fa.clean.nr.bed.fa.nr99 > all.global.fa.clean.nr.bed.fa.nr99.bed

clean-up-gene-models-local:

	#cd /mnt/data; python /mnt/source/gimme/src/utils/get_transcript_seq.py all.local.fa.clean.nr.bed chick.fa > all.local.fa.clean.nr.bed.fa
	#cd /mnt/data; cd-hit-est -T 0 -d 0 -c 0.99 -M 8000 -i all.local.fa.clean.nr.bed.fa -o all.local.fa.clean.nr.bed.fa.nr99
	python ../source/gimme/src/utils/cdhit_transcript.py all.local.fa.clean.nr.bed all.local.fa.clean.nr.bed.fa.nr99 > all.local.fa.clean.nr.bed.fa.nr99.bed

clean-up-gene-models-global-local:

	#cd /mnt/data; python /mnt/source/gimme/src/utils/get_transcript_seq.py all.fa.clean.nr.bed chick.fa > all.fa.clean.nr.bed.fa
	#cd /mnt/data; cd-hit-est -T 0 -d 0 -c 0.99 -M 8000 -i all.fa.clean.nr.bed.fa -o all.fa.clean.nr.bed.fa.nr99
	python ../source/gimme/src/utils/cdhit_transcript.py all.fa.clean.nr.bed all.fa.clean.nr.bed.fa.nr99 > all.fa.clean.nr.bed.fa.nr99.bed

find-unique-transcripts:

	cd /mnt/data; python ~/protocol/assembly-diff.py line6u_local.fa.clean.nr line6u_global.fa.clean.nr
	cd /mnt/data; python ~/protocol/assembly-diff.py line6i_local.fa.clean.nr line6i_global.fa.clean.nr
	cd /mnt/data; python ~/protocol/assembly-diff.py line7u_local.fa.clean.nr line7u_global.fa.clean.nr
	cd /mnt/data; python ~/protocol/assembly-diff.py line7i_local.fa.clean.nr line7i_global.fa.clean.nr
	cd /mnt/data; python ~/protocol/assembly-diff.py line6u_global.fa.clean.nr line6u_local.fa.clean.nr
	cd /mnt/data; python ~/protocol/assembly-diff.py line6i_global.fa.clean.nr line6i_local.fa.clean.nr
	cd /mnt/data; python ~/protocol/assembly-diff.py line7u_global.fa.clean.nr line7u_local.fa.clean.nr
	cd /mnt/data; python ~/protocol/assembly-diff.py line7i_global.fa.clean.nr line7i_local.fa.clean.nr

run-blastx:

	cd /mnt/data; \
	export BLASTDB=/mnt/data/data/blastdb; \
	for input in *.uniq.long; do \
		blastx -evalue 1e-20 -outfmt 5 -query $$input -db mouse.proteins -out $$input.xml; \
	done

## Mouse data
#######################
quality-trim-mouse:

	perl /mnt/source/condetri_v2.1.pl -fastq1=SRR203276_1.fastq -fastq2=SRR203276_2.fastq -sc=33 -prefix=SRR203276 -minlen=50 -cutfirst 10

run-tophat-mouse:

	#cd /mnt/data; bowtie2-build mm9.fa mm9
	cd /mnt/data; tophat2 -p 2 -r 150 -o tophat_mouse mm9 SRR203276_trim1.fastq SRR203276_trim2.fastq

extract-reads-mouse:

	cd tophat_mouse_strand; \
	samtools index accepted_hits.bam; \
	for chr in $$(cat ../protocol/mouse.list.txt); do \
		printf "\textracting %s...\n" $$chr; \
		samtools view -b -o local/$$chr.bam accepted_hits.bam $$chr; \
	done

local-velveth-mouse:

	cd tophat_mouse_strand/local; \
	for chr in *bam; do \
		qsub -v input=$$chr ../../protocol/local_velveth.sh; \
	done

local-velvetg-mouse:

	cd tophat_mouse_strand/local; \
	for chr in *bam; do \
		qsub -v input=$$chr ../../protocol/local_velvetg.sh; \
	done

## mRNA and EST alignment
########################

run-blat-mrna-est:

	#cd /mnt/data; blat -noHead -out=psl -mask=lower -extendThroughN chick_3.2bit mrna.fa mrna.fa.psl
	#cd /mnt/data; blat -noHead -out=psl -mask=lower -extendThroughN chick_3.2bit est.fa est.fa.psl
	cd /mnt/data; sort -k 10 est.fa.psl > est.fa.psl.sorted
	cd /mnt/data; ../source/pslReps -nohead -singleHit est.fa.psl.sorted est.psl.best info
	cd /mnt/data; sort -k 10 mrna.fa.psl > mrna.fa.psl.sorted
	cd /mnt/data; ../source/pslReps -nohead -singleHit mrna.fa.psl.sorted mrna.psl.best info

rsem-prepare-reference:

	cat all.fa.clean.nr.bed.fa | python protocol/fasta-to-gene-list.py > all.fa.clean.nr.bed.txt
	cat all.global.fa.clean.nr.bed.fa | python protocol/fasta-to-gene-list.py > all.global.fa.clean.nr.bed.txt
	cat all.local.fa.clean.nr.bed.fa | python protocol/fasta-to-gene-list.py > all.local.fa.clean.nr.bed.txt
	qsub -v list="all.fa.clean.nr.bed.txt",input="all.fa.clean.nr.bed.fa",sample="all.fa.clean.nr" protocol/rsem_prepare_reference.sh
	qsub -v list="all.global.fa.clean.nr.bed.txt",input="all.global.fa.clean.nr.bed.fa",sample="all.global.fa.clean.nr" protocol/rsem_prepare_reference.sh
	qsub -v list="all.local.fa.clean.nr.bed.txt",input="all.local.fa.clean.nr.bed.fa",sample="all.local.fa.clean.nr" protocol/rsem_prepare_reference.sh

rsem-calculate-expr:

	qsub -v input_read="line6u.se.fq",sample_name="line6u-rsem-all",index="all.fa.clean.nr" \
		protocol/rsem_calculate_expr_single.sh
	qsub -v input_read="line6i.se.fq",sample_name="line6i-rsem-all",index="all.fa.clean.nr" \
		protocol/rsem_calculate_expr_single.sh
	qsub -v input_read="line7u.se.fq",sample_name="line7u-rsem-all",index="all.fa.clean.nr" \
		protocol/rsem_calculate_expr_single.sh
	qsub -v input_read="line7i.se.fq",sample_name="line7i-rsem-all",index="all.fa.clean.nr" \
		protocol/rsem_calculate_expr_single.sh

	qsub -v input_read="line6u.se.fq",sample_name="line6u-rsem-global",index="all.global.fa.clean.nr" \
		protocol/rsem_calculate_expr_single.sh
	qsub -v input_read="line6i.se.fq",sample_name="line6i-rsem-global",index="all.global.fa.clean.nr" \
		protocol/rsem_calculate_expr_single.sh
	qsub -v input_read="line7u.se.fq",sample_name="line7u-rsem-global",index="all.global.fa.clean.nr" \
		protocol/rsem_calculate_expr_single.sh
	qsub -v input_read="line7i.se.fq",sample_name="line7i-rsem-global",index="all.global.fa.clean.nr" \
		protocol/rsem_calculate_expr_single.sh

	qsub -v input_read="line6u.se.fq",sample_name="line6u-rsem-local",index="all.local.fa.clean.nr" \
		protocol/rsem_calculate_expr_single.sh
	qsub -v input_read="line6i.se.fq",sample_name="line6i-rsem-local",index="all.local.fa.clean.nr" \
		protocol/rsem_calculate_expr_single.sh
	qsub -v input_read="line7u.se.fq",sample_name="line7u-rsem-local",index="all.local.fa.clean.nr" \
		protocol/rsem_calculate_expr_single.sh
	qsub -v input_read="line7i.se.fq",sample_name="line7i-rsem-local",index="all.local.fa.clean.nr" \
		protocol/rsem_calculate_expr_single.sh

cufflinks:

	/mnt/source/cufflinks2/cufflinks -o line6u_cufflinks line6u_tophat/accepted_hits.bam
	/mnt/source/cufflinks2/cufflinks -o line6i_cufflinks line6u_tophat/accepted_hits.bam
	/mnt/source/cufflinks2/cufflinks -o line7u_cufflinks line6u_tophat/accepted_hits.bam
	/mnt/source/cufflinks2/cufflinks -o line7i_cufflinks line6u_tophat/accepted_hits.bam
	ls line*cufflinks/transcripts.gtf >> cufflinks.list
	cuffmerge -o cuffmerge -s chick.fa -p 2 cufflinks.list

gimme-assembly-cufflinks:

	python /mnt/source/gimme/src/utils/gff2bed.py cuffmerge/transcripts.gtf > cuffmerge/transcripts.bed
	python /mnt/source/gimme/src/gimme.py all.fa.clean.nr.bed cuffmerge/transcripts.bed > all.fa.clean.nr.cuff.bed

gimme-assembly-cufflinks-mrna:

	python /mnt/source/gimme/src/gimme.py mrna.psl.best  all.fa.clean.nr.bed cuffmerge/transcripts.bed > all.fa.clean.nr.cuff.mrna.bed

filter-low-isopct:

	python protocol/filter-low-isopct.py 1.0 all.fa.clean.nr.bed *rsem-all.isoforms.results > all.fa.clean.nr.flt.bed
	python protocol/filter-low-isopct.py 1.0 all.global.fa.clean.nr.bed *rsem-global.isoforms.results > all.global.fa.clean.nr.flt.bed
	python protocol/filter-low-isopct.py 1.0 all.local.fa.clean.nr.bed *rsem-local.isoforms.results > all.local.fa.clean.nr.flt.bed
	python protocol/filter-low-isopct.py 1.0 all.fa.clean.nr.cuff.bed *rsem-all-cuff.isoforms.results > all.fa.clean.nr.cuff.flt.bed
	python protocol/filter-low-isopct.py 1.0 all.fa.clean.nr.cuff.mrna.bed *rsem-all-cuff-mrna.isoforms.results > all.fa.clean.nr.cuff.mrna.flt.bed

rsem-prepare-reference-cuff-and-asm:

	#python ~/gimme/src/utils/get_transcript_seq.py all.fa.clean.nr.cuff.bed chick.fa > all.fa.clean.nr.cuff.bed.fa
	#cat all.fa.clean.nr.cuff.bed.fa | python protocol/fasta-to-gene-list.py > all.fa.clean.nr.cuff.bed.txt
	qsub -v list="all.fa.clean.nr.cuff.bed.txt",input="all.fa.clean.nr.cuff.bed.fa",sample="all.fa.clean.nr.cuff" protocol/rsem_prepare_reference.sh

rsem-prepare-reference-cuff-and-asm-and-mrna:

	cat all.fa.clean.nr.cuff.mrna.bed.fa | python protocol/fasta-to-gene-list.py > all.fa.clean.nr.cuff.mrna.bed.txt
	qsub -v list="all.fa.clean.nr.cuff.mrna.bed.txt",input="all.fa.clean.nr.cuff.mrna.bed.fa",sample="all.fa.clean.nr.cuff.mrna" protocol/rsem_prepare_reference.sh

rsem-prepare-reference-pasa:

	cat pasa_all_4.assemblies.fasta | python ../protocol/fasta-to-gene-list.py > pasa_all_4.assemblies.bed.txt
	qsub -v list="pasa_all_4.assemblies.bed.txt",input="pasa_all_4.assemblies.fasta",sample="pasa-assemblies" ../protocol/rsem_prepare_reference.sh

rsem-calculate-expr-2:

	qsub -v input_read="line6u.se.fq",sample_name="line6u-rsem-all-cuff",index="all.fa.clean.nr.cuff" \
		protocol/rsem_calculate_expr_single.sh
	qsub -v input_read="line6i.se.fq",sample_name="line6i-rsem-all-cuff",index="all.fa.clean.nr.cuff" \
		protocol/rsem_calculate_expr_single.sh
	qsub -v input_read="line7u.se.fq",sample_name="line7u-rsem-all-cuff",index="all.fa.clean.nr.cuff" \
		protocol/rsem_calculate_expr_single.sh
	qsub -v input_read="line7i.se.fq",sample_name="line7i-rsem-all-cuff",index="all.fa.clean.nr.cuff" \
		protocol/rsem_calculate_expr_single.sh

	qsub -v input_read="line6u.se.fq",sample_name="line6u-rsem-all-cuff-mrna",index="all.fa.clean.nr.cuff.mrna" \
		protocol/rsem_calculate_expr_single.sh
	qsub -v input_read="line6i.se.fq",sample_name="line6i-rsem-all-cuff-mrna",index="all.fa.clean.nr.cuff.mrna" \
		protocol/rsem_calculate_expr_single.sh
	qsub -v input_read="line7u.se.fq",sample_name="line7u-rsem-all-cuff-mrna",index="all.fa.clean.nr.cuff.mrna" \
		protocol/rsem_calculate_expr_single.sh
	qsub -v input_read="line7i.se.fq",sample_name="line7i-rsem-all-cuff-mrna",index="all.fa.clean.nr.cuff.mrna" \
		protocol/rsem_calculate_expr_single.sh
