quality-trim:

	qsub -v "fastq1=reads/line6i.se.fq" $(protocol)/quality_trim_se.sh
	qsub -v "fastq1=reads/line6u.se.fq" $(protocol)/quality_trim_se.sh
	qsub -v "fastq1=reads/line7u.se.fq" $(protocol)/quality_trim_se.sh
	qsub -v "fastq1=reads/line7i.se.fq" $(protocol)/quality_trim_se.sh

velveth-global:

	qsub -v "outdir=line6u_global,input=reads/line6u.fq_trim.fastq" \
		$(protocol)/velveth_job.sh
	qsub -v "outdir=line6i_global,input=reads/line6i.fq_trim.fastq" \
		$(protocol)/velveth_job.sh
	qsub -v "outdir=line7u_global,input=reads/line7u.fq_trim.fastq" \
		$(protocol)/velveth_job.sh
	qsub -v "outdir=line7i_global,input=reads/line7i.fq_trim.fastq" \
		$(protocol)/velveth_job.sh

velvetg-global:

	qsub -v "inputdir=line6u_global" $(protocol)/velvetg_job.sh
	qsub -v "inputdir=line6i_global" $(protocol)/velvetg_job.sh
	qsub -v "inputdir=line7u_global" $(protocol)/velvetg_job.sh
	qsub -v "inputdir=line7i_global" $(protocol)/velvetg_job.sh

oases-global:

	qsub -v "inputdir=line6u_global" $(protocol)/oases_job.sh
	qsub -v "inputdir=line6i_global" $(protocol)/oases_job.sh
	qsub -v "inputdir=line7u_global" $(protocol)/oases_job.sh
	qsub -v "inputdir=line7i_global" $(protocol)/oases_job.sh

velveth-M:

	qsub $(protocol)/velveth_M.sh

velvetg-M:

	qsub $(protocol)/velvetg_M.sh

oases-M:

	qsub $(protocol)/oases_M.sh

run-blat-global-two-kmers:

	# cat line*global_21/transcripts.fa > global_k21.fa
	# cat line*global_31/transcripts.fa > global_k31.fa
	# seqclean global_k21.fa -c 8
	# seqclean global_k31.fa -c 8
	# cd-hit-est -T 8 -d 0 -c 1.0 -M 8000 -i global_k21.fa.clean -o global_k21.fa.clean.nr
	# cd-hit-est -T 8 -d 0 -c 1.0 -M 8000 -i global_k31.fa.clean -o global_k31.fa.clean.nr

	qsub -v "input=global_k21.fa.clean.nr,\
		index=chick_3.2bit" $(protocol)/blat_job.sh

	qsub -v "input=global_k31.fa.clean.nr,\
		index=chick_3.2bit" $(protocol)/blat_job.sh

run-blat-oases-M:

	cd global_merged; seqclean transcripts.fa -c 8

	cd global_merged; \
		qsub -v "input=transcripts.fa.clean,\
		index=../chick_3.2bit" $(protocol)/blat_job.sh

tophat:

	qsub -v "input=reads/line6u.fq_trim.fastq,outdir=line6u_tophat,\
		index=gga3" $(protocol)/tophat.sh
	qsub -v "input=reads/line6i.fq_trim.fastq,outdir=line6i_tophat,\
		index=gga3" $(protocol)/tophat.sh
	qsub -v "input=reads/line7u.fq_trim.fastq,outdir=line7u_tophat,\
		index=gga3" $(protocol)/tophat.sh
	qsub -v "input=reads/line7i.fq_trim.fastq,outdir=line7i_tophat,\
		index=gga3" $(protocol)/tophat.sh

# extract-reads:
# 
# 	cd line6u_tophat; samtools index accepted_hits.bam
# 	cd line6i_tophat; samtools index accepted_hits.bam
# 	cd line7u_tophat; samtools index accepted_hits.bam
# 	cd line7i_tophat; samtools index accepted_hits.bam
# 
# 	cd line6u_tophat; \
# 		for chr in $$(cat ../chromosomes.list); do \
# 			printf "\textracting %s...\n" $$chr; \
# 			samtools view -b -o $$chr.bam accepted_hits.bam $$chr; done
# 
# 	cd line6i_tophat; \
# 		for chr in $$(cat ../chromosomes.list); do \
# 			printf "\textracting %s...\n" $$chr; \
# 			samtools view -b -o $$chr.bam accepted_hits.bam $$chr; done
# 
# 	cd line7u_tophat; \
# 		for chr in $$(cat ../chromosomes.list); do \
# 			printf "\textracting %s...\n" $$chr; \
# 			samtools view -b -o $$chr.bam accepted_hits.bam $$chr; done
# 
# 	cd line7i_tophat; \
# 		for chr in $$(cat ../chromosomes.list); do \
# 			printf "\textracting %s...\n" $$chr; \
# 			samtools view -b -o $$chr.bam accepted_hits.bam $$chr; done

local-velveth:

	cd line6u_tophat; \
		qsub -v "outdir=assembly,input=accepted_hits.bam" \
			$(protocol)/local_velveth.sh
	cd line6i_tophat; \
		qsub -v "outdir=assembly,input=accepted_hits.bam" \
			$(protocol)/local_velveth.sh
	cd line7u_tophat; \
		qsub -v "outdir=assembly,input=accepted_hits.bam" \
			$(protocol)/local_velveth.sh
	cd line7i_tophat; \
		qsub -v "outdir=assembly,input=accepted_hits.bam" \
			$(protocol)/local_velveth.sh


local-velvetg:

	cd line6u_tophat; \
		qsub -v "outdir=assembly" $(protocol)/velvetg_job.sh
	cd line6i_tophat; \
		qsub -v "outdir=assembly" $(protocol)/velvetg_job.sh
	cd line7u_tophat; \
		qsub -v "outdir=assembly" $(protocol)/velvetg_job.sh
	cd line7i_tophat; \
		qsub -v "outdir=assembly" $(protocol)/velvetg_job.sh

local-oases:

	cd line6u_tophat; \
		qsub -v "inputdir=assembly" $(protocol)/oases_job.sh
	cd line6i_tophat; \
		qsub -v "inputdir=assembly" $(protocol)/oases_job.sh
	cd line7u_tophat; \
		qsub -v "inputdir=assembly" $(protocol)/oases_job.sh
	cd line7i_tophat; \
		qsub -v "inputdir=assembly" $(protocol)/oases_job.sh

combine-local-assembly-transcripts:

	cd line6u_tophat; \
	for d in assembly_??; do \
	if [ -e $$d/transcripts.fa ]; then \
		python $(gimmedir)/src/utils/rename_fasta.py \
			$$d/transcripts.fa line6u_local_$$d >> ../line6u_local.fa; \
	fi; done

	cd line6i_tophat; \
	for d in assembly_??; do \
	if [ -e $$d/transcripts.fa ]; then \
		python $(gimmedir)/src/utils/rename_fasta.py \
		$$d/transcripts.fa line6i_local_$$d >> ../line6i_local.fa; \
	fi; done

	cd line7u_tophat; \
	for d in assembly_??; do \
	if [ -e $$d/transcripts.fa ]; then \
		python $(gimmedir)/src/utils/rename_fasta.py \
		$$d/transcripts.fa line7u_local_$$d >> ../line7u_local.fa; \
	fi; done

	cd line7i_tophat; \
	for d in assembly_??; do \
	if [ -e $$d/transcripts.fa ]; then \
		python $(gimmedir)/src/utils/rename_fasta.py \
		$$d/transcripts.fa line7i_local_$$d >> ../line7i_local.fa; \
	fi; done

combine-global-assembly-transcripts:

	for d in line6u_global_*; do \
		python $(gimmedir)/src/utils/rename_fasta.py \
		$$d/transcripts.fa $$d >> line6u_global.fa; \
	done

	for d in line6i_global_*; do \
		python $(gimmedir)/src/utils/rename_fasta.py \
		$$d/transcripts.fa $$d >> line6i_global.fa; \
	done

	for d in line7u_global_*; do \
		python $(gimmedir)/src/utils/rename_fasta.py \
		$$d/transcripts.fa $$d >> line7u_global.fa; \
	done

	for d in line7i_global_*; do \
		python $(gimmedir)/src/utils/rename_fasta.py \
		$$d/transcripts.fa $$d >> line7i_global.fa; \
	done

clean-transcripts:

	seqclean line6u_local.fa -c 16
	seqclean line6i_local.fa -c 16
	seqclean line7u_local.fa -c 16
	seqclean line7i_local.fa -c 16
	seqclean line6u_global.fa -c 16
	seqclean line6i_global.fa -c 16
	seqclean line7u_global.fa -c 16
	seqclean line7i_global.fa -c 16

remove-redundant-transcripts:

	cd-hit-est -T 8 -d 0 -c 1.0 -M 8000 -i line6u_local.fa.clean -o line6u_local.fa.clean.nr
	cd-hit-est -T 8 -d 0 -c 1.0 -M 8000 -i line6i_local.fa.clean -o line6i_local.fa.clean.nr
	cd-hit-est -T 8 -d 0 -c 1.0 -M 8000 -i line7u_local.fa.clean -o line7u_local.fa.clean.nr
	cd-hit-est -T 8 -d 0 -c 1.0 -M 8000 -i line7i_local.fa.clean -o line7i_local.fa.clean.nr

	cd-hit-est -T 8 -d 0 -c 1.0 -M 8000 -i line6u_global.fa.clean -o line6u_global.fa.clean.nr
	cd-hit-est -T 8 -d 0 -c 1.0 -M 8000 -i line6i_global.fa.clean -o line6i_global.fa.clean.nr
	cd-hit-est -T 8 -d 0 -c 1.0 -M 8000 -i line7u_global.fa.clean -o line7u_global.fa.clean.nr
	cd-hit-est -T 8 -d 0 -c 1.0 -M 8000 -i line7i_global.fa.clean -o line7i_global.fa.clean.nr

remove-redundant-all-assembly:

	cat line*_global.fa line*_local.fa > all_assembly.fa.clean
	qsub -v "c=1.0,input=all_assembly.fa.clean,output=all_assembly.clean.nr" \
		$(protocol)/cdhit_job.sh

remove-redundant-local-assembly:

	cat line*_local.fa > local_assembly.fa.clean
	qsub -v "c=1.0,input=local_assembly.fa.clean,output=local_assembly.clean.nr" \
		$(protocol)/cdhit_job.sh

remove-redundant-global-assembly:

	cat line*_global.fa > global_assembly.fa.clean
	qsub -v "c=1.0,input=global_assembly.fa.clean,output=global_assembly.clean.nr" \
		$(protocol)/cdhit_job.sh

run-blat-global-assembly:

	qsub -v "input=global_assembly.clean.nr,index=chick_3.2bit" $(protocol)/blat_job.sh

run-blat-local-assembly:

	qsub -v "input=local_assembly.clean.nr,index=chick_3.2bit" $(protocol)/blat_job.sh

run-blat-all-assembly:

	qsub -v "input=all_assembly.clean.nr,index=chick_3.2bit" $(protocol)/blat_job.sh

construct-gene-models-global:

	sort -k 10 global_assembly.clean.nr.psl > global_assembly.clean.nr.psl.sorted
	pslReps -nohead -singleHit global_assembly.clean.nr.psl.sorted \
		global_assembly.clean.nr.psl.best info

	qsub -v input="global_assembly.clean.nr.psl.best,\
		output=global_assembly_models.bed,gimme_dir=$(gimmedir)/src/" \
		$(protocol)/run_gimme.sh

construct-gene-models-local:

	sort -k 10 local_assembly.clean.nr.psl > local_assembly.clean.nr.psl.sorted
	pslReps -nohead -singleHit local_assembly.clean.nr.psl.sorted \
		local_assembly.clean.nr.psl.best info

	qsub -v input="local_assembly.clean.nr.psl.best,\
		output=local_assembly_models.bed,gimme_dir=$(gimmedir)/src/" \
		$(protocol)/run_gimme.sh

construct-gene-models-global-local:

	sort -k 10 all_assembly.clean.nr.psl > all_assembly.clean.nr.psl.sorted
	pslReps -nohead -singleHit all_assembly.clean.nr.psl.sorted \
		all_assembly.clean.nr.psl.best info

	qsub -v input="all_assembly.clean.nr.psl.best,\
		output=all_assembly_models.bed,gimme_dir=$(gimmedir)/src/" \
		$(protocol)/run_gimme.sh

# clean-up-gene-models-global-local:

	#cd /mnt/data; python /mnt/source/gimme/src/utils/get_transcript_seq.py all.fa.clean.nr.bed chick.fa > all.fa.clean.nr.bed.fa
	#cd /mnt/data; cd-hit-est -T 0 -d 0 -c 0.99 -M 8000 -i all.fa.clean.nr.bed.fa -o all.fa.clean.nr.bed.fa.nr99
	# python ../source/gimme/src/utils/cdhit_transcript.py all.fa.clean.nr.bed all.fa.clean.nr.bed.fa.nr99 > all.fa.clean.nr.bed.fa.nr99.bed

find-unique-transcripts:

	# python $(protocol)/assembly-diff.py line6u_local.fa.clean.nr line6u_global.fa.clean.nr
	python $(protocol)/assembly-diff.py line6i_local.fa.clean.nr line6i_global.fa.clean.nr
	python $(protocol)/assembly-diff.py line7u_local.fa.clean.nr line7u_global.fa.clean.nr
	python $(protocol)/assembly-diff.py line7i_local.fa.clean.nr line7i_global.fa.clean.nr
	python $(protocol)/assembly-diff.py line6u_global.fa.clean.nr line6u_local.fa.clean.nr
	python $(protocol)/assembly-diff.py line6i_global.fa.clean.nr line6i_local.fa.clean.nr
	python $(protocol)/assembly-diff.py line7u_global.fa.clean.nr line7u_local.fa.clean.nr
	python $(protocol)/assembly-diff.py line7i_global.fa.clean.nr line7i_local.fa.clean.nr

run-blastx:

	cd /mnt/data; \
	export BLASTDB=/mnt/data/data/blastdb; \
	for input in *.uniq.long; do \
		blastx -evalue 1e-20 -outfmt 5 -query $$input -db mouse.proteins -out $$input.xml; \
	done

rsem-prepare-reference-global-local:

	python $(gimmedir)/src/utils/get_transcript_seq.py \
		all_assembly_models.bed chick.fa | grep -v DEBUG >  all_assembly_models.bed.fa
	cat all_assembly_models.bed.fa | python $(protocol)/fasta-to-gene-list.py \
		> all_assembly_models.bed.fa.txt

	qsub -v list="all_assembly_models.bed.fa.txt,\
		input=all_assembly_models.bed.fa,sample=all-assembly-models" \
		$(protocol)/rsem_prepare_reference.sh

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

	qsub -v "outdir=line6u,input=line6u_tophat/accepted_hits.bam" $(protocol)/cufflinks.sh
	qsub -v "outdir=line6i,input=line6i_tophat/accepted_hits.bam" $(protocol)/cufflinks.sh
	qsub -v "outdir=line7u,input=line7u_tophat/accepted_hits.bam" $(protocol)/cufflinks.sh
	qsub -v "outdir=line7i,input=line7i_tophat/accepted_hits.bam" $(protocol)/cufflinks.sh

cufflinks-merge:

	ls line*cuff/transcripts.gtf >> cufflinks.list
	cuffmerge -o cuffmerge -s chick.fa -p 8 cufflinks.list

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

## mRNA and EST alignment

run-blat-mrna-est:

	qsub -v "index=chick_3.2bit,input=mrna.fa" $(protocol)/blat_job.sh
	qsub -v "index=chick_3.2bit,input=est.fa" $(protocol)/blat_job.sh

sort-mrna-est-alignments:

	sort -k 10 est.fa.psl > est.fa.psl.sorted
	pslReps -nohead -singleHit est.fa.psl.sorted est.psl.best info
	sort -k 10 mrna.fa.psl > mrna.fa.psl.sorted
	pslReps -nohead -singleHit mrna.fa.psl.sorted mrna.psl.best info

## Mouse data

quality-trim-mouse:

	perl /mnt/source/condetri_v2.1.pl -fastq1=SRR203276_1.fastq -fastq2=SRR203276_2.fastq -sc=33 -prefix=SRR203276 -minlen=50 -cutfirst 10

interleave-mouse-reads:

	shuffleSequences_fastq.pl SRR203276_trim1.fastq SRR203276_trim2.fastq mouse-paired.fastq
	
velveth-mouse:

	qsub -v "pe_input=mouse-paired.fastq,\
		se_input=SRR203276_trim_unpaired.fastq" $(protocol)/velveth_mouse.sh

velvetg-mouse:

	qsub $(protocol)/velvetg_mouse.sh

oases-mouse:

	qsub $(protocol)/oases_mouse.sh

tophat-mouse:

	qsub -v "outdir=tophat_mouse,index=mm9,\
		left=SRR203276_trim1.fastq,right=SRR203276_trim2.fastq,\
		unpaired=SRR203276_trim_unpaired.fastq" \
		$(protocol)/tophat_mouse.sh

run-tophat-mouse-old:

	qsub -v "outdir=tophat_mouse,index=mm9,\
		left=SRR203276_trim1.fastq,right=SRR203276_trim2.fastq,\
		unpaired=SRR203276_trim_unpaired.fastq" tophat_mouse.sh

cufflinks-mouse:

	qsub -v "outdir=mouse_cufflinks,input=tophat_mouse/accepted_hits.bam" $(protocol)/cufflinks.sh

# extract-reads-mouse:
# 
# 	cd tophat_mouse; \
# 	samtools index accepted_hits.bam; \
# 	for chr in $$(cat $(protocol)/mouse.list.txt); do \
# 		printf "\textracting %s...\n" $$chr; \
# 		samtools view -b -o $$chr.bam accepted_hits.bam $$chr; \
# 	done

local-velveth-mouse:

	cd tophat_mouse; \
		qsub -v "outdir=assembly,input=accepted_hits.bam" \
		$(protocol)/local_velveth_mouse.sh

local-velvetg-mouse:

	cd tophat_mouse; \
		qsub -v outdir=assembly $(protocol)/local_velvetg_mouse.sh

local-oases-mouse:

	cd tophat_mouse; \
			qsub -v dir=assembly $(protocol)/oases_mouse.sh

combine-global-mouse-assembly-transcripts:

	for d in mouse_global_*; do \
		python $(gimmedir)/src/utils/rename_fasta.py \
		$$d/transcripts.fa $$d >> mouse_global.fa; \
	done

combine-local-mouse-assembly-transcripts:

	cd tophat_mouse; \
		for d in assembly_??; do \
			python $(gimmedir)/src/utils/rename_fasta.py \
			$$d/transcripts.fa mouse_local_$$d >> ../mouse_local.fa; \
		done

clean-remove-redundant-mouse-transcripts:

	seqclean mouse_local.fa -c 8
	seqclean mouse_global.fa -c 8

	cat mouse_local.fa.clean mouse_global.fa.clean > mouse_assembly.fa.clean
	qsub -v "c=1.0,input=mouse_assembly.fa.clean,output=mouse_assembly.clean.nr" \
		$(protocol)/cdhit_job.sh

run-blat-all-mouse-assembly:

	qsub -v "input=mouse_assembly.clean.nr,index=mm9.fa" $(protocol)/blat_job.sh

construct-gene-models-mouse:

	# sort -k 10 mouse_assembly.clean.nr.psl > mouse_assembly.clean.nr.psl.sorted
	# pslReps -nohead -singleHit mouse_assembly.clean.nr.psl.sorted \
	# 	mouse_assembly.clean.nr.psl.best info

	# grep -w $$chr mouse_assembly.clean.nr.psl.best > mouse_$$chr.psl; \

	for chr in $$(cat $(protocol)/mouse.list.txt); do \
		qsub -v input="mouse_$$chr.psl,output=mouse_$$chr.bed,\
		gimme_dir=$(gimmedir)/src/" $(protocol)/run_gimme.sh; \
	done
