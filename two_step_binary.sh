#!/bin/bash

#PBS -N tsb
#PBS -l walltime=00:05:00
#PBS -t 0-10

cd $PBS_O_WORKDIR

module load python
module load numpy

python two_step_binary.py $PBS_ARRAYID
