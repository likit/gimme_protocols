#!/bin/sh -login
#PBS -l nodes=1:ppn=9,mem=12gb,walltime=24:00:00
#PBS -M preeyano@msu.edu
#PBS -m abe
#PBS -N BLAST+

module load BLAST+
cd ${PBS_O_WORKDIR}

${prog} -num_threads 8 -outfmt 5 -db ${db} -out ${out} -query ${query}
