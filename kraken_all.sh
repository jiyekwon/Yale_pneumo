#!/bin/bash
#SBATCH --partition scavenge
#SBATCH --job-name=kraken_mini
#SBATCH --error=kraken_mini.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=kraken_%j.log

# Path to your NEW mini database folder
module load miniconda
eval "$(conda shell.bash hook)"
conda activate /vast/palmer/scratch/weinberger_daniel/jk2666/zephyr_pneumo_project/envs/pneumo_capsule


BASE_DIR="/vast/palmer/scratch/weinberger_daniel/jk2666/zephyr_pneumo_project"
DB_PATH="${BASE_DIR}/databases"
DATA_DIR="${BASE_DIR}/data"
SRR_FILE="${BASE_DIR}/srr_list.txt"

ACCESSION=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $SRR_FILE)
echo "Processing $ACCESSION (Task ID: $SLURM_ARRAY_TASK_ID)"

# download
fasterq-dump $ACCESSION --threads 8 --outdir $DATA_DIR
fastp -i ${DATA_DIR}/${ACCESSION}.fastq -o ${DATA_DIR}/${ACCESSION}_trimmed.fastq --thread 8

# 5. Run Kraken2
kraken2 --db "$DB_PATH" \
        --threads 8 \
        --report "${DATA_DIR}/${ACCESSION}_report.txt" \
        "${DATA_DIR}/${ACCESSION}_trimmed.fastq" > "${DATA_DIR}/${ACCESSION}.kraken"