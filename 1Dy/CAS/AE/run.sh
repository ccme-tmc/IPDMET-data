#!/bin/bash
#SBATCH -p amd
#SBATCH --qos=low
#SBATCH -J lizw
#SBATCH --output=jobid%j-%N.out
#SBATCH --error=jobid%j-%N.out
#SBATCH --ntasks-per-node=16

python -u main.py
