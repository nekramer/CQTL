#!/bin/bash
#SBATCH -J PosToRsids
#SBATCH -t 12:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -p general

module load python/3.9.6
python3 scripts/posToRsids.py

