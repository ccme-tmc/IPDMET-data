#!/bin/bash
#SBATCH -p cp6
#SBATCH -o testall.out
#SBATCH -J pyscf-test
#SBATCH -n 56

#Task
python CePhMe.py -i ./


