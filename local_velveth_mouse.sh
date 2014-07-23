#!/bin/sh -login
#PBS -l nodes=1:ppn=1,mem=48gb,walltime=24:00:00
#PBS -m abe
#PBS -N Velveth_local_mouse
#PBS -M preeyano@msu.edu

cd ${PBS_O_WORKDIR}

velveth ${dir} 21,33,2 -bam -shortPaired -strand_specific ${input}
