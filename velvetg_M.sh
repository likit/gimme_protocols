#!/bin/sh -login
#PBS -l nodes=1:ppn=1,mem=64gb,walltime=8:00:00
#PBS -M preeyano@msu.edu
#PBS -m abe
#PBS -A ged-intel11
#PBS -N Velvetg_global_merge_gimme

cd ${PBS_O_WORKDIR}
velvetg global_merged -read_trkg yes -conserveLong yes
