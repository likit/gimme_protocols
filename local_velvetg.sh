#!/bin/sh -login
#PBS -l nodes=1:ppn=1,mem=48gb,walltime=24:00:00
#PBS -m abe
#PBS -N Velvetg_local_mouse
#PBS -t 21,23,25,27,29,31

cd ${PBS_O_WORKDIR}

dir=$(basename ${input} .bam)
/mnt/home/preeyano/velvet_1.2.03/velvetg ${dir}_asm_${PBS_ARRAYID} -read_trkg yes -ins_length 300
