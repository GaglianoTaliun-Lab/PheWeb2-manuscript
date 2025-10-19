#!/bin/bash
#SBATCH --account=def-gsarah
#SBATCH --job-name=gwas_combined
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --array=0-451
#SBATCH --time=00:30:00

source .venv/bin/activate

GWAS_DIR="/lustre06/project/6060121/CLSA_PheWeb_shared/PheWeb_v2_final_data/generated-by-pheweb/parsed"  
SCRIPT="get_significant_loci.py"          
OUTPUT_DIR="loci_combined"
MIN_AF=0.01

FILES=($(ls ${GWAS_DIR}/*.combined* | grep -v '\.interaction-'))
GWAS_FILE="${FILES[$SLURM_ARRAY_TASK_ID]}"
#GWAS_FILE="/home/justb11/scratch/Thesis_work/copy_generated_by_pheweb/generated-by-pheweb/parsed/BLD_MPV_COM.all.combined"

# Run script
python $SCRIPT \
    -g "$GWAS_FILE" \
    -af "$MIN_AF" \
    -o "$OUTPUT_DIR"
