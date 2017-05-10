#!/bin/bash

#PBS -N tsbc
#PBS -l walltime=00:30:00
#PBS -t 1-100

cd $PBS_O_WORKDIR

module load python/3.5.1

python two_step_binary_cont.py
