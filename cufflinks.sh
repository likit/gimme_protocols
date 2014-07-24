#!/bin/sh -login
#PBS -l nodes=1:ppn=5,mem=24gb,walltime=24:00:00
#PBS -M preeyano@msu.edu
#PBS -m abe
#PBS -N Cufflinks_gimme${PBS_JOBID}

cd ${PBS_O_WORKDIR}
module load cufflinks/2.0.0
cufflinks -p 4 -o ${outdir}_cuff ${input}
