#!/bin/sh -login
#PBS -l nodes=1:ppn=1,mem=48gb,walltime=24:00:00
#PBS -M preeyano@msu.edu
#PBS -m abe
#PBS -N Gimme_${PBS_JOBID}

module load bxPython
module load pygr
module load matplotlib

cd ${PBS_O_WORKDIR}
python ${gimme_dir}/gimme.py ${input} > ${output} 2>${PBS_JOBID}-gimme.log
