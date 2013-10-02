#!/bin/sh -login
#PBS -l nodes=1:ppn=1,mem=48gb,walltime=24:00:00
#PBS -m abe
#PBS -N Velveth_line6u_global

cd /mnt/ls12/preeyanon/gimme/
velveth line6u_global 21,33,2 -fastq -short 6147JAAXX_2_1_pf_trim.fastq
