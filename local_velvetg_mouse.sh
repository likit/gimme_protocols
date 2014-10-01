#!/bin/sh -login
#PBS -l nodes=1:ppn=1,mem=48gb,walltime=24:00:00
#PBS -m abe
#PBS -N Velvetg_local_mouse

cd ${PBS_O_WORKDIR}

velvetg ${outdir} -read_trkg yes -ins_length 300
