#!/bin/sh -login
#PBS -l nodes=1:ppn=8,mem=24gb,walltime=24:00:00
#PBS -M preeyano@msu.edu
#PBS -m abe
#PBS -N Tophat_single_${PBS_JOBID}

module load bowtie/1.0.0
module load TopHat/1.3.1
cd ${PBS_O_WORKDIR}

tophat -p 7 -o ${outdir} ${index} ${input}
