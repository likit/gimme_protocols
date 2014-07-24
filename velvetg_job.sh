#!/bin/sh -login
#PBS -l nodes=1:ppn=1,mem=64gb,walltime=24:00:00
#PBS -M preeyano@msu.edu
#PBS -m abe
#PBS -A ged-intel11
#PBS -N Velvetg_${PBS_JOBID}
#PBS -t 21,23,25,27,29,31

cd ${PBS_O_WORKDIR}
velvetg ${outdir}_${PBS_ARRAYID} -read_trkg yes
