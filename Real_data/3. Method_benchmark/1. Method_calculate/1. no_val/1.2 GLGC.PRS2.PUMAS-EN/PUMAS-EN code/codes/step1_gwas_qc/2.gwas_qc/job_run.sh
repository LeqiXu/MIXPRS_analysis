#!/bin/bash
#SBATCH --mem=8G
#SBATCH -p scavenge,day,week
#SBATCH -t 1:00:00
#SBATCH -c 1

# --- Load Modules ---
module load R/4.3.0-foss-2020b

# --- Paths ---
CODE_DIR="/home/yd357/pi_paths/pi_zhao/GWAS_Subsample_JointPRS/v2_PUMAS-EN/PUMAS-EN/codes/step1_gwas_qc/2.gwas_qc"
GWAS_DIR="/home/yd357/pi_paths/pi_zhao/GWAS_Subsample_JointPRS/data/summary_data/X_Wing"
FRQ_DIR="/home/yd357/pi_paths/pi_zhao/GWAS_Subsample_JointPRS/v2_PUMAS-EN/PUMAS-EN/data/LD_ref/A1A2.match/2_make_LD_genotype"
OUTPUT_DIR="/home/yd357/pi_paths/pi_zhao/GWAS_Subsample_JointPRS/v2_PUMAS-EN/PUMAS-EN/data/GWAS"
mkdir -p "${OUTPUT_DIR}"

# --- Read Parameters ---
IFS=',' read -r TRAIT POP N_SAMPLES GROUP < <(sed -n "${SLURM_ARRAY_TASK_ID}p" ${CODE_DIR}/job_params.txt)

# --- Define File Paths ---
INPUT_FILE="${GWAS_DIR}/${TRAIT}_${POP}_inter.txt"
FRQ_FILE="${FRQ_DIR}/${GROUP}/${POP}/qced/1KG.LD500.PRSCS.1kg.A1A2match.freq.frq"

# --- Run QC ---
Rscript /home/yd357/pi_paths/pi_zhao/GWAS_Subsample_JointPRS/Project21_PUMAS-ensemble/tool/code/gwas_qc.R \
  --file_path "${INPUT_FILE}" \
  --frq_path "${FRQ_FILE}" \
  --output_path "${OUTPUT_DIR}" \
  --snp SNP \
  --a1 A1 \
  --a2 A2 \
  --stat BETA \
  --p P \
  --chr CHR \
  --bp BP \
  --n.total "${N_SAMPLES}"
