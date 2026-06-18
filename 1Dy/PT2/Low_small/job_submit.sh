#!/bin/bash
#SBATCH -p cp4
#SBATCH -o testall.out
#SBATCH -J pyscf-test
#SBATCH -n 56

#Task
python DyCp2.py -i ./


