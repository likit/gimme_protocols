#!/bin/sh -login
#PBS -l nodes=1:ppn=5,mem=24gb,walltime=24:00:00
#PBS -M preeyano@msu.edu
#PBS -m abe
#PBS -N Tophat_paired_${PBS_JOBID}

module load bowtie2/2.1.0
cd ${PBS_O_WORKDIR}

~/tophat-2.0.9.Linux_x86_64/tophat -r 150 -p 4 -o ${outdir} ${index} ${left} ${right},${unpaired}
