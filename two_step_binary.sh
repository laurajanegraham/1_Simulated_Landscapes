#!/bin/bash

#PBS -N tsb
#PBS -l walltime=00:15:00
#PBS -t 1-100

cd $PBS_O_WORKDIR

module load python
module load numpy

python two_step_binary.py
