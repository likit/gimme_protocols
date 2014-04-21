#!/bin/sh -login
#PBS -l nodes=1:ppn=1,mem=48gb,walltime=24:00:00
#PBS -m abe
#PBS -N Velveth_local_mouse

cd ${PBS_O_WORKDIR}

/mnt/home/preeyano/velvet_1.2.03/velveth $(basename ${input} .bam)_asm 21,33,2 -bam -shortPaired -strand_specific ${input}
