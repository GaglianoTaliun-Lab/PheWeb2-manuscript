#!/bin/bash
#SBATCH --account=def-gsarah
#SBATCH --job-name=loci_compare
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=20G
#SBATCH --time=02:50:00

source .venv/bin/activate

python compare_sig_loci.py \
    -m total.all.male_loci_list.txt \
    -f total.all.female_loci_list.txt \
    -c total.all.combined_loci_list.txt