#!/bin/sh -login
#PBS -l nodes=1:ppn=1,mem=64gb,walltime=48:00:00
#PBS -m abe
#PBS -v DATAPATH
#PBS -N Oases_line7i_global
#PBS -t 21,23,25,27,29,31

cd $DATAPATH
oases line7i_global_${PBS_ARRAYID} -unused_reads yes
