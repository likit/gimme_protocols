#!/bin/sh -login
#PBS -l nodes=1:ppn=1,mem=16gb,walltime=12:00:00
#PBS -m abe
#PBS -N condetri_line6i
#PBS -A ged-intel11
#PBS -v DATAPATH

cd $DATAPATH
condetri_v2.1.pl -fastq1=6147JAAXX_3_1_pf.fastq -sc=33 -cutfirst 10
