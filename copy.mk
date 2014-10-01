copypath=preeyano@hpc.msu.edu:/mnt/ls12/preeyanon/gimme-paper

all:

	# scp -C $(copypath)/*uniq results/
	# scp -C $(copypath)/global*k??*all_sp.txt results/
	# scp -C $(copypath)/line*global*clean.nr results/
	# scp -C $(copypath)/line*local*clean.nr results/
	# scp -C $(copypath)/global_merged/transcripts.fa.clean.psl.best_all_sp.txt \
	# 	results/oases-M.psl.best.all_sp.txt
	# scp -C $(copypath)/Gallus*all_sp.txt results/
	# scp -C $(copypath)/cuffmerge/merged.bed_all_sp.txt \
	# 	results/cufflinks.bed_all_sp.txt
	# scp -C $(copypath)/cuffmerge-ref/merged.bed_all_sp.txt \
	# 	results/cufflinks-ref.bed_all_sp.txt
	# scp -C $(copypath)/cuffmerge-ref/merged.bed \
	# 	results/cufflinks-ref.bed
	# scp -C $(copypath)/cuffmerge/merged.bed \
	# 	results/cufflinks.bed
	# scp -C $(copypath)/*bed results/

	# scp -C $(copypath)/*.out results/
	# scp -C $(copypath)/*long results/
	# scp -C $(copypath)/*gimme*counts results/
	# scp -C $(copypath)/cuffmerge/*counts results/
	# scp -C $(copypath)/global*all_sp.txt results/
	# scp -C $(copypath)/local*all_sp.txt results/
	# scp -C $(copypath)/all_assembly_models.bed.faa results/
	# scp -C $(copypath)/mouse_assembly_models.bed results/
	# scp -C $(copypath)/mouse_assembly_models.bed_all_sp.txt results/
	# scp -C $(copypath)/mouse_cuff/transcripts.select.bed_all_sp.txt \
	# 	results/mouse_cufflinks.bed_all_sp.txt
	# scp -C $(copypath)/Mus_musculus.NCBIM37.64.bed_all_sp.txt results/
	# scp -C $(copypath)/line6u_untrim/line6u.profile results/
	# scp -C $(copypath)/line6i_untrim/line6i.profile results/
	# scp -C $(copypath)/line7u_untrim/line7u.profile results/
	# scp -C $(copypath)/line7i_untrim/line7i.profile results/
	#
	# scp -C $(copypath)/cuffmerge/merged.bed.faa.blastp.out \
	# 	results/cufflinks.bed.faa.blastp.out
	# scp -C $(copypath)/cuffmerge/merged.bed.faa results/cufflinks.bed.faa
	# scp -C $(copypath)/all_assembly_cufflinks_models.bed.faa.blastp.out \
	# 	results/all_assembly_cufflinks_models.bed.faa.blastp.out
	scp -C $(copypath)/all_assembly_cufflinks_models.bed.faa \
		results/all_assembly_cufflinks_models.bed.faa
