#!/bin/sh -login
#PBS -l nodes=1:ppn=1,mem=16gb,walltime=12:00:00
#PBS -m abe
#PBS -N condetri_line7i
#PBS -A ged-intel11

cd /mnt/ls12/preeyanon/gimme/
condetri_v2.1.pl -fastq1=6147JAAXX_7_1_pf.fastq -sc=33 -cutfirst 10
