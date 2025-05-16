#!/bin/bash
#SBATCH --mem=16G
#SBATCH --partition=scavenge,day,week,pi_zhao
#SBATCH -t 1:30:00
#SBATCH -c 1

module load R/4.3.0-foss-2020b

# --- Configuration ---
BASE_PROJECT_DIR="/home/yd357/pi_paths/pi_zhao/GWAS_Subsample_JointPRS/v3/concatenate"
CODE_DIR_STEP4="${BASE_PROJECT_DIR}/codes/array_job/step4_filtered_pumas_eval"
JOB_PARAMS_FILE="${BASE_PROJECT_DIR}/codes/array_job/step4_pumas_eval/job_params_step4_pumas_eval.txt"
FILTERED_LIST_DIR="${CODE_DIR_STEP4}/filtered_model_lists"

PUMAS_SCRIPT_DIR="/home/yd357/pi_paths/pi_zhao/GWAS_Subsample_JointPRS/Project21_PUMAS-ensemble/tool/code"
PUMAS_EVAL_SCRIPT_NAME="PUMAS-ensemble.evaluation.R"
LD_REF_BASE_DIR="/home/yd357/pi_paths/pi_zhao/GWAS_Subsample_JointPRS/v2_PUMAS-EN/PUMAS-EN/data/LD_ref/A1A2.match/2_make_LD_genotype/GLGC"
SUBSAMPLING_DATA_BASE_DIR="/home/yd357/pi_paths/pi_zhao/GWAS_Subsample_JointPRS/v2_PUMAS-EN/PUMAS-EN/data/subsampling/PUMA-ensemble.subsampling"
SUBSAMPLING_ENSEMBLE_METHOD_FOLDER="EN"

SEPARATED_BASE_DIR="/home/yd357/pi_paths/pi_zhao/GWAS_Subsample_JointPRS/v3/concatenate/results/separate"
SUB_BETA_DIR_BASE="${SEPARATED_BASE_DIR}/sub_beta"
FULL_BETA_DIR_BASE="${SEPARATED_BASE_DIR}/full_beta"
PUMAS_EVAL_OUTPUT_BASE_DIR="${BASE_PROJECT_DIR}/results/pumas_evaluation_output_filtered"
K_FOLDS=4

# --- Read Job Parameters ---
IFS=',' read -r TARGET_TRAIT TARGET_POPULATION < <(sed -n "${SLURM_ARRAY_TASK_ID}p" "${JOB_PARAMS_FILE}")

PUMAS_SCRIPT_TRAIT_NAME="${TARGET_TRAIT}_${TARGET_POPULATION}_inter"
LD_REF_PATH="${LD_REF_BASE_DIR}/${TARGET_POPULATION}/qced/1KG.LD500.PRSCS.1kg.A1A2match"
XTY_STATS_PATH="${SUBSAMPLING_DATA_BASE_DIR}/${SUBSAMPLING_ENSEMBLE_METHOD_FOLDER}/${PUMAS_SCRIPT_TRAIT_NAME}/"
FINAL_PUMAS_OUTPUT_DIR="${PUMAS_EVAL_OUTPUT_BASE_DIR}/${TARGET_TRAIT}/${TARGET_POPULATION}/"
SUB_BETA_DIR="${SUB_BETA_DIR_BASE}/${PUMAS_SCRIPT_TRAIT_NAME}/"
FULL_BETA_DIR="${FULL_BETA_DIR_BASE}/${PUMAS_SCRIPT_TRAIT_NAME}/"
MODEL_LIST_FILTERED="${FILTERED_LIST_DIR}/${TARGET_TRAIT}_${TARGET_POPULATION}.filtered_models.txt"

mkdir -p "${FINAL_PUMAS_OUTPUT_DIR}"

# --- Get Filtered PRS Method List ---
PRS_METHODS_STR=$(paste -sd, "${MODEL_LIST_FILTERED}")
NUM_FILTERED_MODELS=$(wc -l < "${MODEL_LIST_FILTERED}")

# --- Run PUMAS Evaluation ---
cd "${PUMAS_SCRIPT_DIR}"
Rscript "${PUMAS_EVAL_SCRIPT_NAME}" \
    --k "${K_FOLDS}" \
    --ensemble "EN" \
    --ref_path "${LD_REF_PATH}" \
    --trait_name "${PUMAS_SCRIPT_TRAIT_NAME}" \
    --prs_method "${PRS_METHODS_STR}" \
    --xty_path "${XTY_STATS_PATH}" \
    --stats_path "${XTY_STATS_PATH}" \
    --weight_path "${SUB_BETA_DIR}" \
    --full_weight_path "${FULL_BETA_DIR}" \
    --output_path "${FINAL_PUMAS_OUTPUT_DIR}/" \
    --threads 1

# --- Output Validation ---
OUTPUT_WEIGHT_FILE="${FINAL_PUMAS_OUTPUT_DIR}/${PUMAS_SCRIPT_TRAIT_NAME}.ensemble.weights.txt"
if [ -s "${OUTPUT_WEIGHT_FILE}" ]; then
    NUM_OUTPUT_COLS=$(head -n 1 "${OUTPUT_WEIGHT_FILE}" | awk -F'\t' '{print NF}')
    EXPECTED_COLS=$((NUM_FILTERED_MODELS + 3))
fi
