#!/bin/bash

#PBS -N tsb
#PBS -l walltime=00:10:00
#PBS -t 1-100

cd $PBS_O_WORKDIR

module load python/3.5.1

python farmland_birds.py
