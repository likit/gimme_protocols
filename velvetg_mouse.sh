#!/bin/sh -login
#PBS -l nodes=1:ppn=1,mem=128gb,walltime=24:00:00
#PBS -M preeyano@msu.edu
#PBS -m abe
#PBS -A ged-intel11
#PBS -N Velvetg_global_${PBS_JOBID}

cd ${PBS_O_WORKDIR}
velvetg mouse_global_27 -read_trkg yes -ins_length 300
