#!/bin/bash
#SBATCH --mem=24G
#SBATCH -p scavenge,day,week 
#SBATCH -t 1:00:00
#SBATCH -c 1

# --- Load Required Modules ---
module load R/4.3.0-foss-2020b

# --- Set Paths ---
CODE_DIR="/home/yd357/pi_paths/pi_zhao/GWAS_Subsample_JointPRS/v2_PUMAS-EN/PUMAS-EN/codes/step2_subsampling"
GWAS_DIR="/home/yd357/pi_paths/pi_zhao/GWAS_Subsample_JointPRS/v2_PUMAS-EN/PUMAS-EN/data/GWAS"
LD_DIR="/home/yd357/pi_paths/pi_zhao/GWAS_Subsample_JointPRS/v2_PUMAS-EN/PUMAS-EN/data/LD_ref/RData"
SCRIPT_DIR="/home/yd357/pi_paths/pi_zhao/GWAS_Subsample_JointPRS/Project21_PUMAS-ensemble/tool/code"

# --- Read Parameters ---
IFS=',' read -r TRAIT POP ENSEMBLE_METHOD < <(sed -n "${SLURM_ARRAY_TASK_ID}p" ${CODE_DIR}/job_params.txt)

# --- Define Paths ---
TRAIT_NAME="${TRAIT}_${POP}_inter"
LD_FILE="${LD_DIR}/1000G_${POP}.RData"
BASE_OUTPUT_DIR="/home/yd357/pi_paths/pi_zhao/GWAS_Subsample_JointPRS/v2_PUMAS-EN/PUMAS-EN/data/subsampling/PUMA-ensemble.subsampling/${ENSEMBLE_METHOD}"
OUTPUT_DIR="${BASE_OUTPUT_DIR}/${TRAIT_NAME}"
mkdir -p "${OUTPUT_DIR}"

# --- Run Subsampling ---
Rscript "${SCRIPT_DIR}/PUMA-ensemble.subsampling.R" \
  --k 4 \
  --partitions 0.6,0.2,0.1,0.1 \
  --trait_name "${TRAIT_NAME}" \
  --ensemble "${ENSEMBLE_METHOD}" \
  --gwas_path "${GWAS_DIR}" \
  --ld_path "${LD_FILE}" \
  --output_path "${OUTPUT_DIR}/" \
  --parallel \
  --threads 1
