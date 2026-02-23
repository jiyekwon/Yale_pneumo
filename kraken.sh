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
DB_PATH="/vast/palmer/scratch/weinberger_daniel/jk2666/zephyr_pneumo_project/databases"

# Run Kraken2 on your trimmed single-end file
kraken2 --db $DB_PATH \
        --threads 8 \
        --report ./data/SRR36951856_mini_report.txt \
        ./data/SRR36951856_trimmed.fastq > ./data/SRR36951856.kraken