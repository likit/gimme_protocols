#!/bin/sh -login
#PBS -l nodes=1:ppn=1,mem=48gb,walltime=24:00:00
#PBS -M preeyano@msu.edu
#PBS -m abe
#PBS -N Oases_gimme${PBS_JOBID}
#PBS -t 21,23,25,27,29,31

cd ${PBS_O_WORKDIR}

oases ${inputdir}_${PBS_ARRAYID}
