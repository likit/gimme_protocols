#!/bin/sh -login
#PBS -l nodes=1:ppn=4,mem=24gb,walltime=12:00:00
#PBS -M preeyano@msu.edu
#PBS -m abe
#PBS -N Bowtie

module load bowtie/1.0.0
cd ${PBS_O_WORKDIR}
bowtie -M 100 -S -t -p 4 ${index} ${input} ${output}
