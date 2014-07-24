#!/bin/sh -login
#PBS -l nodes=1:ppn=1,mem=64gb,walltime=8:00:00
#PBS -m abe
#PBS -N Velveth_global_merge_gimme

cd ${PBS_O_WORKDIR}
velveth global_merged 27 -long line*global_*/transcripts.fa
