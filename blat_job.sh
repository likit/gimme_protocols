#!/bin/sh -login
#PBS -l nodes=1:ppn=1,mem=24gb,walltime=24:00:00
#PBS -M preeyano@msu.edu
#PBS -m abe
#PBS -N Blat_${PBS_JOBID}

cd ${PBS_O_WORKDIR}
blat -noHead -extendThroughN -mask=lower -out=psl ${index} ${input} ${input}.psl
