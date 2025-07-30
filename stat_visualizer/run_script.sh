#!/bin/bash
#SBATCH --account=def-gsarah
#SBATCH --time=06:00:00
#SBATCH --mem-per-cpu=10G
#SBATCH --cpus-per-task=10
#SBATCH --ntasks=1
#SBATCH --nodes=1

module load StdEnv/2020 python/3.11.5
#source myenv/bin/activate
python3 run_multiple_files.py
