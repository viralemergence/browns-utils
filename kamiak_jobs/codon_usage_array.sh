#!/bin/bash
#SBATCH --partition=remi
#SBATCH --job-name=codon_usage_array
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alexander.brown@wsu.edu
#SBATCH --time=7-00:00:00
#SBATCH --array=0-60:1%10
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G

CDS_DIR="/data/lab/seifert/CDS_domestic/sequences_cds/"
CODON_USAGE_DIR="/data/lab/seifert/codon_usage/"

SCRIPT="/home/alexander.brown/browns-utils/codon_usage_calculator.py"

python3 $SCRIPT \
    --cds $CDS_DIR \
    --outdir $CODON_USAGE_DIR \
    --file_index $SLURM_ARRAY_TASK_ID