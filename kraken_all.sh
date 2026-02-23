#!/bin/bash
#SBATCH --partition scavenge
#SBATCH --job-name=kraken_pneumo
#SBATCH --array=1-6
#SBATCH --error=kraken_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=kraken_%A_%a.log

module load miniconda
eval "$(conda shell.bash hook)"
conda activate /nfs/roberts/scratch/pi_dmw63/jk2666/zephyr_pneumo_project/envs/pneumo_capsule


BASE_DIR="/nfs/roberts/scratch/pi_dmw63/jk2666/zephyr_pneumo_project"
DB_PATH="${BASE_DIR}/databases"
DATA_DIR="${BASE_DIR}/data"
SRR_FILE="${BASE_DIR}/srr_list.txt"

ACCESSION=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $SRR_FILE)
echo "Processing $ACCESSION (Task ID: $SLURM_ARRAY_TASK_ID)"

# download
fasterq-dump --split-files $ACCESSION --threads 8 --outdir $DATA_DIR
fastp -i ${DATA_DIR}/${ACCESSION}.fastq -o ${DATA_DIR}/${ACCESSION}_trimmed.fastq --thread 8

# 5. Run Kraken2
kraken2 --db "$DB_PATH" \
        --threads 8 \
        --report "${DATA_DIR}/${ACCESSION}_report.txt" \
        "${DATA_DIR}/${ACCESSION}_trimmed.fastq" > "${DATA_DIR}/${ACCESSION}.kraken"