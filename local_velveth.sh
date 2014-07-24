#!/bin/sh -login
#PBS -l nodes=1:ppn=1,mem=64gb,walltime=24:00:00
#PBS -m abe
#PBS -N Velveth_global_gimme
#PBS -M preeyano@msu.edu

cd ${PBS_O_WORKDIR}
velveth ${outdir} 21,33,2 -bam -short ${input}
