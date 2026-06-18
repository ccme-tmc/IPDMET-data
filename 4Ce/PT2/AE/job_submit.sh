#!/bin/bash
#SBATCH -p cp6bm
#SBATCH -o testall.out
#SBATCH -J pyscf-test
#SBATCH -n 56

#Task
python CeNP2O2.py -i ./


