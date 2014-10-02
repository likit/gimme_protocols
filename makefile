quality-trim:

	qsub -v "input=reads/line6i.se.fq" $(protocol)/quality_trim_se_job.sh
	qsub -v "input=reads/line6u.se.fq" $(protocol)/quality_trim_se_job.sh
	qsub -v "input=reads/line7u.se.fq" $(protocol)/quality_trim_se_job.sh
	qsub -v "input=reads/line7i.se.fq" $(protocol)/quality_trim_se_job.sh

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

run-blat-oases-M:

	cd global_merged; seqclean transcripts.fa -c 8

	cd global_merged; \
		qsub -v "input=transcripts.fa.clean,\
		index=../chick_3.2bit" $(protocol)/blat_job.sh

sort-oases-M-alignments:

	cd global_merged; \
		sort -k 10 transcripts.fa.clean.psl > transcripts.fa.clean.psl.sorted; \
		pslReps -nohead -singleHit transcripts.fa.clean.psl.sorted \
		transcripts.fa.clean.psl.best info

tophat:

	qsub -v "input=reads/line6u.fq_trim.fastq,outdir=line6u_tophat,\
		index=gga3" $(protocol)/tophat.sh
	qsub -v "input=reads/line6i.fq_trim.fastq,outdir=line6i_tophat,\
		index=gga3" $(protocol)/tophat.sh
	qsub -v "input=reads/line7u.fq_trim.fastq,outdir=line7u_tophat,\
		index=gga3" $(protocol)/tophat.sh
	qsub -v "input=reads/line7i.fq_trim.fastq,outdir=line7i_tophat,\
		index=gga3" $(protocol)/tophat.sh

extract-reads:

	cd line6u_tophat; samtools index accepted_hits.bam
	cd line6i_tophat; samtools index accepted_hits.bam
	cd line7u_tophat; samtools index accepted_hits.bam
	cd line7i_tophat; samtools index accepted_hits.bam

	cd line6u_tophat; \
		for chr in $$(cat ../chromosomes.list); do \
			printf "\textracting %s...\n" $$chr; \
			samtools view -b -o $$chr.bam accepted_hits.bam $$chr; done

	cd line6i_tophat; \
		for chr in $$(cat ../chromosomes.list); do \
			printf "\textracting %s...\n" $$chr; \
			samtools view -b -o $$chr.bam accepted_hits.bam $$chr; done

	cd line7u_tophat; \
		for chr in $$(cat ../chromosomes.list); do \
			printf "\textracting %s...\n" $$chr; \
			samtools view -b -o $$chr.bam accepted_hits.bam $$chr; done

	cd line7i_tophat; \
		for chr in $$(cat ../chromosomes.list); do \
			printf "\textracting %s...\n" $$chr; \
			samtools view -b -o $$chr.bam accepted_hits.bam $$chr; done

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
	for d in chr*_??; do \
	if [ -e $$d/transcripts.fa ]; then \
		python $(gimmedir)/src/utils/rename_fasta.py \
			$$d/transcripts.fa line6u_local_$$d >> ../line6u_local.fa; \
	fi; done

	cd line6i_tophat; \
	for d in chr*_??; do \
	if [ -e $$d/transcripts.fa ]; then \
		python $(gimmedir)/src/utils/rename_fasta.py \
		$$d/transcripts.fa line6i_local_$$d >> ../line6i_local.fa; \
	fi; done

	cd line7u_tophat; \
	for d in chr*_??; do \
	if [ -e $$d/transcripts.fa ]; then \
		python $(gimmedir)/src/utils/rename_fasta.py \
		$$d/transcripts.fa line7u_local_$$d >> ../line7u_local.fa; \
	fi; done

	cd line7i_tophat; \
	for d in chr*_??; do \
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

	cat line*_global.fa.clean.nr line*_local.fa.clean.nr > all_assembly.fa.clean
	qsub -v "c=1.0,input=all_assembly.fa.clean,output=all_assembly.clean.nr" \
		$(protocol)/cdhit_job.sh

remove-redundant-local-assembly:

	cat line*_local.fa.clean.nr > local_assembly.fa.clean
	qsub -v "c=1.0,input=local_assembly.fa.clean,output=local_assembly.clean.nr" \
		$(protocol)/cdhit_job.sh

remove-redundant-global-assembly:

	cat line*_global.fa.clean.nr > global_assembly.fa.clean
	qsub -v "c=1.0,input=global_assembly.fa.clean,output=global_assembly.clean.nr" \
		$(protocol)/cdhit_job.sh

run-blat-global-assembly:

	for f in global_assembly.clean.nr_*.fa; do \
		qsub -v "input=$$f,index=chick_3.2bit" $(protocol)/blat_job.sh; \
	done

run-blat-local-assembly:

	for f in local_assembly.clean.nr_*.fa; do \
		qsub -v "input=$$f,index=chick_3.2bit" $(protocol)/blat_job.sh; \
	done

run-blat-all-assembly:

	for f in all_assembly.clean.nr_*.fa; do \
		qsub -v "input=$$f,index=chick_3.2bit" $(protocol)/blat_job.sh; \
	done

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

find-unique-transcripts:

	python $(protocol)/assembly-diff.py line6u_local.fa.clean.nr line6u_global.fa.clean.nr
	python $(protocol)/assembly-diff.py line6i_local.fa.clean.nr line6i_global.fa.clean.nr
	python $(protocol)/assembly-diff.py line7u_local.fa.clean.nr line7u_global.fa.clean.nr
	python $(protocol)/assembly-diff.py line7i_local.fa.clean.nr line7i_global.fa.clean.nr
	python $(protocol)/assembly-diff.py line6u_global.fa.clean.nr line6u_local.fa.clean.nr
	python $(protocol)/assembly-diff.py line6i_global.fa.clean.nr line6i_local.fa.clean.nr
	python $(protocol)/assembly-diff.py line7u_global.fa.clean.nr line7u_local.fa.clean.nr
	python $(protocol)/assembly-diff.py line7i_global.fa.clean.nr line7i_local.fa.clean.nr

get-long-unique-regions:

	python $(gimmedir)/src/utils/sizeSelect.py \
		line6u_global.fa.clean.nr.uniq 300 > line6u_global.fa.clean.nr.uniq.long
	python $(gimmedir)/src/utils/sizeSelect.py \
		line6i_global.fa.clean.nr.uniq 300 > line6i_global.fa.clean.nr.uniq.long
	python $(gimmedir)/src/utils/sizeSelect.py \
		line7u_global.fa.clean.nr.uniq 300 > line7u_global.fa.clean.nr.uniq.long
	python $(gimmedir)/src/utils/sizeSelect.py \
		line7i_global.fa.clean.nr.uniq 300 > line7i_global.fa.clean.nr.uniq.long

	python $(gimmedir)/src/utils/sizeSelect.py \
		line6u_local.fa.clean.nr.uniq 300 > line6u_local.fa.clean.nr.uniq.long
	python $(gimmedir)/src/utils/sizeSelect.py \
		line6i_local.fa.clean.nr.uniq 300 > line6i_local.fa.clean.nr.uniq.long
	python $(gimmedir)/src/utils/sizeSelect.py \
		line7u_local.fa.clean.nr.uniq 300 > line7u_local.fa.clean.nr.uniq.long
	python $(gimmedir)/src/utils/sizeSelect.py \
		line7i_local.fa.clean.nr.uniq 300 > line7i_local.fa.clean.nr.uniq.long

cufflinks:

	qsub -v "outdir=line6u,input=line6u_tophat/accepted_hits.bam" $(protocol)/cufflinks.sh
	qsub -v "outdir=line6i,input=line6i_tophat/accepted_hits.bam" $(protocol)/cufflinks.sh
	qsub -v "outdir=line7u,input=line7u_tophat/accepted_hits.bam" $(protocol)/cufflinks.sh
	qsub -v "outdir=line7i,input=line7i_tophat/accepted_hits.bam" $(protocol)/cufflinks.sh

cufflinks-merge:

	ls line*cuff/transcripts.gtf > cufflinks.list
	cuffmerge -o cuffmerge -s chick.fa -p 8 cufflinks.list

# cufflinks-merge-ref:
# 
# 	ls line*cuff/transcripts.gtf > cufflinks.list
# 	cuffmerge -g Gallus_gallus.WASHUC2.64.gtf -o cuffmerge-ref -s chick.fa -p 8 cufflinks.list
# 
construct-gene-models-global-local-cufflinks:

	python $(protocol)/gff2bed.py cuffmerge/merged.gtf > cuffmerge/merged.bed
	qsub -v input="all_assembly_models.bed cuffmerge/merged.bed,\
		output=all_assembly_cufflinks_models.bed,gimme_dir=$(gimmedir)/src/" \
		$(protocol)/run_gimme.sh

construct-gene-models-global-local-cufflinks-ensembl:

	python $(protocol)/gff2bed.py cuffmerge-ref/merged.gtf > cuffmerge-ref/merged.bed
	qsub -v input="all_assembly_models.bed cuffmerge-ref/merged.bed,\
		output=all_assembly_cufflinks_ensembl_models.bed,gimme_dir=$(gimmedir)/src/" \
		$(protocol)/run_gimme.sh

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

	qsub -v "outdir=mouse,input=tophat_mouse/accepted_hits.bam" $(protocol)/cufflinks.sh

extract-reads-mouse:

	cd tophat_mouse; \
	samtools index accepted_hits.bam; \
	for chr in $$(cat $(protocol)/mouse.list.txt); do \
		printf "\textracting %s...\n" $$chr; \
		samtools view -b -o $$chr.bam accepted_hits.bam $$chr; \
	done

local-velveth-mouse:

	cd tophat_mouse; \
		for chr in chr*.bam; do \
			velveth $$(basename $$chr .bam) 27 -bam -shortPaired $$chr -strand_specific; \
		done

local-velvetg-mouse:

	cd tophat_mouse; \
		for chr in $$(cat $(protocol)/mouse.list.txt); do \
			qsub -v outdir=$$chr $(protocol)/local_velvetg_mouse.sh; \
		done

local-oases-mouse:

	cd tophat_mouse; \
		for chr in $$(cat $(protocol)/mouse.list.txt); do \
			oases $$chr -ins_length 300; \
		done

combine-global-mouse-assembly-transcripts:

	for d in mouse_global_27; do \
		python $(gimmedir)/src/utils/rename_fasta.py \
		$$d/transcripts.fa $$d >> mouse_global.fa; \
	done

combine-local-mouse-assembly-transcripts:

	cd tophat_mouse; \
		for chr in $$(cat $(protocol)/mouse.list.txt); do \
			cat $$chr/transcripts.fa >> ../mouse_local.fa; \
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

	sort -k 10 mouse_assembly.clean.nr.psl > mouse_assembly.clean.nr.psl.sorted
	pslReps -nohead -singleHit mouse_assembly.clean.nr.psl.sorted \
		mouse_assembly.clean.nr.psl.best info
	
	# split positive and negative strand transcripts
	
	awk '$$9=="+"' mouse_assembly.clean.nr.psl.best > mouse_assembly.clean.nr.psl.best.pos
	awk '$$9=="-"' mouse_assembly.clean.nr.psl.best > mouse_assembly.clean.nr.psl.best.neg

	qsub -v input="mouse_assembly.clean.nr.psl.best.neg,\
		output=mouse_assembly_models.neg.bed,gimme_dir=$(gimmedir)/src/" \
		$(protocol)/run_gimme_mouse.sh

	qsub -v input="mouse_assembly.clean.nr.psl.best.pos,\
		output=mouse_assembly_models.pos.bed,gimme_dir=$(gimmedir)/src/" \
		$(protocol)/run_gimme_mouse.sh

	cat mouse_assembly_models.pos.bed mouse_assembly_models.neg.bed > mouse_assembly_models.bed
	python $(protocol)/select-mouse-genes.py $(protocol)/mouse.list.txt mouse_assembly_models.bed > tmp
	mv tmp mouse_assembly_models.bed


#############################
# Junctions analysis
############################

run-blat-global-two-kmers:

	cat line*global_21/transcripts.fa > global_k21.fa
	cat line*global_31/transcripts.fa > global_k31.fa
	seqclean global_k21.fa -c 8
	seqclean global_k31.fa -c 8
	cd-hit-est -T 8 -d 0 -c 1.0 -M 8000 -i global_k21.fa.clean -o global_k21.fa.clean.nr
	cd-hit-est -T 8 -d 0 -c 1.0 -M 8000 -i global_k31.fa.clean -o global_k31.fa.clean.nr

	qsub -v "input=global_k21.fa.clean.nr,\
		index=chick_3.2bit" $(protocol)/blat_job.sh

	qsub -v "input=global_k31.fa.clean.nr,\
		index=chick_3.2bit" $(protocol)/blat_job.sh

sort-blat-two-kmers:

	sort -k 10 global_k21.fa.clean.nr.psl > global_k21.fa.clean.nr.psl.sorted
	pslReps -nohead -singleHit global_k21.fa.clean.nr.psl.sorted \
		global_k21.fa.clean.nr.psl.best info

	sort -k 10 global_k31.fa.clean.nr.psl > global_k31.fa.clean.nr.psl.sorted
	pslReps -nohead -singleHit global_k31.fa.clean.nr.psl.sorted \
		global_k31.fa.clean.nr.psl.best info

find-splice-junctions:

	python $(gimmedir)/src/utils/compare_junction.py --all -p global_k21.fa.clean.nr.psl.best
	python $(gimmedir)/src/utils/compare_junction.py --all -p global_k31.fa.clean.nr.psl.best
	cd global_merged; \
		python $(gimmedir)/src/utils/compare_junction.py --all -p \
		transcripts.fa.clean.psl.best
	python $(gimmedir)/src/utils/compare_junction.py --all -p est.psl.best
	python $(gimmedir)/src/utils/compare_junction.py --all -p mrna.psl.best
	python $(gimmedir)/src/utils/compare_junction.py --all -p \
		global_assembly.clean.nr.psl.best
	cd cuffmerge; \
		python $(protocol)/gff2bed.py merged.gtf > merged.bed; \
		python $(gimmedir)/src/utils/compare_junction.py --all -b merged.bed
	cd cuffmerge-ref; \
		python $(protocol)/gff2bed.py merged.gtf > merged.bed; \
		python $(gimmedir)/src/utils/compare_junction.py --all -b merged.bed
	python $(gimmedir)/src/utils/compare_junction.py --all -b all_assembly_models.bed
	python $(gimmedir)/src/utils/compare_junction.py --all -b all_assembly_cufflinks_models.bed
	python $(gimmedir)/src/utils/compare_junction.py --all -b all_assembly_cufflinks_ensembl_models.bed
	
	python $(protocol)/src/utils/gff2bed.py Gallus_gallus.WASHUC2.64.gtf > Gallus_gallus.WASHUC2.64.bed
	python $(gimmedir)/src/utils/compare_junction.py --all -b Mus_musculus.NCBIM37.64.bed
	python $(gimmedir)/src/utils/compare_junction.py --all -b Gallus_gallus.WASHUC2.64.bed
	
	python $(gimmedir)/src/utils/compare_junction.py --all -b mouse_assembly_models.bed
	python $(protocol)/gff2bed.py mouse_cuff/transcripts.gtf > mouse_cuff/transcripts.bed
	cd mouse_cuff; \
		python $(protocol)/select-mouse-genes.py $(protocol)/mouse.list.txt \
		transcripts.bed > transcripts.select.bed; \
		python $(gimmedir)/src/utils/compare_junction.py --all -b transcripts.select.bed

#####################
# Read mappings
####################

ensembl-map:

	qsub -v "input=reads/line6u.fq_trim.fastq,output=line6u.ensembl.sam,\
		index=Gallus_gallus.WASHUC2.64.cdna.all.fa" \
		$(protocol)/bowtie_job.sh
	qsub -v "input=reads/line6i.fq_trim.fastq,output=line6i.ensembl.sam,\
		index=Gallus_gallus.WASHUC2.64.cdna.all.fa" \
		$(protocol)/bowtie_job.sh
	qsub -v "input=reads/line7u.fq_trim.fastq,output=line7u.ensembl.sam,\
		index=Gallus_gallus.WASHUC2.64.cdna.all.fa" \
		$(protocol)/bowtie_job.sh
	qsub -v "input=reads/line7i.fq_trim.fastq,output=line7i.ensembl.sam,\
		index=Gallus_gallus.WASHUC2.64.cdna.all.fa" \
		$(protocol)/bowtie_job.sh

	qsub -v "left=reads/line6u.1_trim1.fastq,\
		right=reads/line6u.1_trim2.fastq,output=line6u.ensembl.pe.sam,\
		index=Gallus_gallus.WASHUC2.64.cdna.all.fa" \
		$(protocol)/bowtie_pe_job.sh

	qsub -v "left=reads/line6i.1_trim1.fastq,\
		right=reads/line6i.1_trim2.fastq,output=line6i.ensembl.pe.sam,\
		index=Gallus_gallus.WASHUC2.64.cdna.all.fa" $(protocol)/bowtie_pe_job.sh

	qsub -v "left=reads/line7u.1_trim1.fastq,\
		right=reads/line7u.1_trim2.fastq,output=line7u.ensembl.pe.sam,\
		index=Gallus_gallus.WASHUC2.64.cdna.all.fa" $(protocol)/bowtie_pe_job.sh

	qsub -v "left=reads/line7i.1_trim1.fastq,\
		right=reads/line7i.1_trim2.fastq,output=line7i.ensembl.pe.sam,\
		index=Gallus_gallus.WASHUC2.64.cdna.all.fa" $(protocol)/bowtie_pe_job.sh

build-bowtie-index-gimme-models:

	python $(gimmedir)/src/utils/get_transcript_seq.py \
		all_assembly_models.bed chick.fa > all_assembly_models.bed.fa
	bowtie-build all_assembly_models.bed.fa all_assembly_models.bed.fa

build-bowtie-index-cufflinks-models:

	python $(gimmedir)/src/utils/get_transcript_seq.py \
		cuffmerge/merged.bed chick.fa > cuffmerge/merged.bed.fa
	cd cuffmerge; \
		bowtie-build merged.bed.fa merged.bed.fa

cufflinks-map:

	cd cuffmerge; \
	qsub -v "input=../reads/line6u.fq_trim.fastq,output=line6u.cufflinks.sam,\
		index=merged.bed.fa" \
		$(protocol)/bowtie_job.sh
	cd cuffmerge; \
	qsub -v "input=../reads/line6i.fq_trim.fastq,output=line6i.cufflinks.sam,\
		index=merged.bed.fa" \
		$(protocol)/bowtie_job.sh
	cd cuffmerge; \
	qsub -v "input=../reads/line7u.fq_trim.fastq,output=line7u.cufflinks.sam,\
		index=merged.bed.fa" \
		$(protocol)/bowtie_job.sh
	cd cuffmerge; \
	qsub -v "input=../reads/line7i.fq_trim.fastq,output=line7i.cufflinks.sam,\
		index=merged.bed.fa" \
		$(protocol)/bowtie_job.sh

	cd cuffmerge; \
	qsub -v "left=../reads/line6u.1_trim1.fastq,\
		right=../reads/line6u.1_trim2.fastq,output=line6u.cufflinks.pe.sam,\
		index=merged.bed.fa" \
		$(protocol)/bowtie_pe_job.sh

	cd cuffmerge; \
	qsub -v "left=../reads/line6i.1_trim1.fastq,\
		right=../reads/line6i.1_trim2.fastq,output=line6i.cufflinks.pe.sam,\
		index=merged.bed.fa" \
		$(protocol)/bowtie_pe_job.sh

	cd cuffmerge; \
	qsub -v "left=../reads/line7u.1_trim1.fastq,\
		right=../reads/line7u.1_trim2.fastq,output=line7u.cufflinks.pe.sam,\
		index=merged.bed.fa" \
		$(protocol)/bowtie_pe_job.sh

	cd cuffmerge; \
	qsub -v "left=../reads/line7i.1_trim1.fastq,\
		right=../reads/line7i.1_trim2.fastq,output=line7i.cufflinks.pe.sam,\
		index=merged.bed.fa" \
		$(protocol)/bowtie_pe_job.sh

translate-cufflinks:

	cd cuffmerge; estscan -M ../gallus.hm -t merged.bed.faa \
		merged.bed.fa > merged.bed.fna

blast-cufflinks-vs-mouse:

	cd cuffmerge; \
	qsub -v "prog=blastp,db=../mouse.protein.faa,\
		query=merged.bed.faa,\
		out=merged.bed.faa.blastp" \
		$(protocol)/blast_job.sh

parse-cufflinks-vs-mouse-blast:

	cd cuffmerge; \
	python $(gimmedir)/src/utils/blast_hits.py \
		merged.bed.faa.blastp 1e-20 > \
		merged.bed.faa.blastp.out

get-sequences-cufflinks-gimme-models:

	python $(gimmedir)/src/utils/get_transcript_seq.py \
		all_assembly_cufflinks_models.bed chick.fa > \
		all_assembly_cufflinks_models.bed.fa

translate-cufflinks-gimme:

	estscan -M gallus.hm -t all_assembly_cufflinks_models.bed.faa \
		all_assembly_cufflinks_models.bed.fa > \
		all_assembly_cufflinks_models.bed.fna

blast-cufflinks-gimme-vs-mouse:

	qsub -v "prog=blastp,db=mouse.protein.faa,\
		query=all_assembly_cufflinks_models.bed.faa,\
		out=all_assembly_cufflinks_models.bed.faa.blastp" \
		$(protocol)/blast_job.sh

parse-cufflinks-gimme-vs-mouse-blast:

	python $(gimmedir)/src/utils/blast_hits.py \
		all_assembly_cufflinks_models.bed.faa.blastp 1e-20 > \
		all_assembly_cufflinks_models.bed.faa.blastp.out

gimme-models-map:

	qsub -v "input=reads/line6u.fq_trim.fastq,output=line6u.gimme.sam,\
		index=all_assembly_models.bed.fa" \
		$(protocol)/bowtie_job.sh
	qsub -v "input=reads/line6i.fq_trim.fastq,output=line6i.gimme.sam,\
		index=all_assembly_models.bed.fa" \
		$(protocol)/bowtie_job.sh
	qsub -v "input=reads/line7u.fq_trim.fastq,output=line7u.gimme.sam,\
		index=all_assembly_models.bed.fa" \
		$(protocol)/bowtie_job.sh
	qsub -v "input=reads/line7i.fq_trim.fastq,output=line7i.gimme.sam,\
		index=all_assembly_models.bed.fa" \
		$(protocol)/bowtie_job.sh

	qsub -v "left=reads/line6u.1_trim1.fastq,\
		right=reads/line6u.1_trim2.fastq,output=line6u.gimme.pe.sam,\
		index=all_assembly_models.bed.fa" \
		$(protocol)/bowtie_pe_job.sh

	qsub -v "left=reads/line6i.1_trim1.fastq,\
		right=reads/line6i.1_trim2.fastq,output=line6i.gimme.pe.sam,\
		index=all_assembly_models.bed.fa" \
		$(protocol)/bowtie_pe_job.sh

	qsub -v "left=reads/line7u.1_trim1.fastq,\
		right=reads/line7u.1_trim2.fastq,output=line7u.gimme.pe.sam,\
		index=all_assembly_models.bed.fa" \
		$(protocol)/bowtie_pe_job.sh

	qsub -v "left=reads/line7i.1_trim1.fastq,\
		right=reads/line7i.1_trim2.fastq,output=line7i.gimme.pe.sam,\
		index=all_assembly_models.bed.fa" \
		$(protocol)/bowtie_pe_job.sh

count-spliced-reads-ensembl:

	for f in *ensembl*sam; do \
		samtools view -b -S $$f -o $$(basename $$f .sam).bam; \
	done
	
	for f in *ensembl*bam; do \
		samtools sort $$f $$(basename $$f .bam).sorted; \
	done
	
	for f in *ensembl*sorted.bam; do samtools index $$f; done
	
	for f in *ensembl*sorted.bam; do \
		python $(gimmedir)/src/utils/count_spliced_reads.py \
		Gallus_gallus.WASHUC2.64.bed $$f > $$f.reads_counts; \
	done

count-spliced-reads-gimme:

	for f in *gimme*sam; do \
		samtools view -b -S $$f -o $$(basename $$f .sam).bam; \
	done

	for f in *gimme*bam; do \
		samtools sort $$f $$(basename $$f .bam).sorted; \
	done

	for f in *gimme*sorted.bam; do samtools index $$f; done

	for f in *gimme*sorted.bam; do \
		python $(gimmedir)/src/utils/count_spliced_reads.py \
		all_assembly_models.bed $$f > $$f.reads_counts; \
	done

count-spliced-reads-cufflinks:

	cd cuffmerge; \
	for f in *sam; do \
		samtools view -b -S $$f -o $$(basename $$f .sam).bam; \
	done

	cd cuffmerge; \
	for f in *bam; do \
		samtools sort $$f $$(basename $$f .bam).sorted; \
	done
	
	cd cuffmerge; \
		for f in *sorted.bam; do samtools index $$f; done
	
	cd cuffmerge; \
		for f in *sorted.bam; do \
			python $(gimmedir)/src/utils/count_spliced_reads.py \
			merged.bed $$f > $$f.reads_counts; \
		done

#################
# Homologs analysis
################

create-mouse-db:

	wget ftp://ftp.ncbi.nlm.nih.gov/refseq/M_musculus/mRNA_Prot/mouse.protein.faa.gz
	gunzip mouse.protein.faa.gz
	makeblastdb -in mouse.protein.faa -dbtype prot -parse_seqids -out mouse.protein.faa
gimme-vs-mouse-blastp:

	# translate cDNA sequences
	estscan -M gallus.hm -t all_assembly_models.bed.faa \
		all_assembly_models.bed.fa > all_assembly_models.bed.fna

	qsub -v "prog=blastp,db=mouse.protein.faa,\
		query=all_assembly_models.bed.faa,\
		out=all_assembly_models.bed.faa.blastp" \
		$(protocol)/blast_job.sh

	python $(gimmedir)/src/utils/blast_hits.py \
		all_assembly_models.bed.faa.blastp 1e-20 > \
		all_assembly_models.bed.faa.blastp.out

run-blastx-unique-regions:

	for input in *.uniq.long; do \
		qsub -v prog=blastx,query=$$input,db=mouse.protein.faa,out=$$input.xml \
		$(protocol)/blast_job.sh; \
	done

find-mouse-match:

	for f in *uniq.long.xml; do \
		python $(protocol)/find_match.py $$f; \
	done

reads-error-profile:

	qsub -v "input=reads/line6u.se.fq.gz,outdir=line6u_untrim,\
		index=gga3" $(protocol)/tophat.sh
	qsub -v "input=reads/line6i.se.fq.gz,outdir=line6i_untrim,\
		index=gga3" $(protocol)/tophat.sh
	qsub -v "input=reads/line7u.se.fq.gz,outdir=line7u_untrim,\
		index=gga3" $(protocol)/tophat.sh
	qsub -v "input=reads/line7i.se.fq.gz,outdir=line7i_untrim,\
		index=gga3" $(protocol)/tophat.sh

	cd line6u_untrim; samtools calmd -b accepted_hits.bam \
		../chick.fa > line6u.bam; \
		samtools index line6u.bam; \
		python $(gimmedir)/src/utils/read_error_profile.py line6u.bam \
		> line6u.profile

	cd line6i_untrim; samtools calmd -b accepted_hits.bam \
		../chick.fa > line6i.bam; \
		samtools index line6i.bam; \
		python $(gimmedir)/src/utils/read_error_profile.py line6i.bam \
		> line6i.profile

	cd line7u_untrim; samtools calmd -b accepted_hits.bam \
		../chick.fa > line7u.bam; \
		samtools index line7u.bam; \
		python $(gimmedir)/src/utils/read_error_profile.py line7u.bam \
		> line7u.profile

	cd line7i_untrim; samtools calmd -b accepted_hits.bam \
		../chick.fa > line7i.bam; \
		samtools index line7i.bam; \
		python $(gimmedir)/src/utils/read_error_profile.py line7i.bam \
		> line7i.profile
