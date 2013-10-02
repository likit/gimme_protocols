#!/bin/sh -login
#PBS -l nodes=1:ppn=1,mem=64gb,walltime=48:00:00
#PBS -m abe
#PBS -v DATAPATH
#PBS -N Velvetg_line6i_global
#PBS -t 21,23,25,27,29,31

cd $DATAPATH
velvetg line6i_global_${PBS_ARRAYID} -read_trkg yes -unused_reads yes
