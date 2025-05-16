#!/bin/bash
#SBATCH --mem=1G
#SBATCH --partition=scavenge,day,week,pi_zhao
#SBATCH -t 0:20:00
#SBATCH -c 1

# --- Configuration ---
BASE_PROJECT_DIR="/home/yd357/pi_paths/pi_zhao/GWAS_Subsample_JointPRS/v3/concatenate"

CODE_DIR_STEP4="${BASE_PROJECT_DIR}/codes/array_job/step4_filtered_pumas_eval"
WORKER_SCRIPT="${CODE_DIR_STEP4}/job_run_filtered_pumas_eval.sh"
WORKER_LOG_DIR="${CODE_DIR_STEP4}/logs/worker_logs"
JOB_PARAMS_FILE="${BASE_PROJECT_DIR}/codes/array_job/step4_pumas_eval/job_params_step4_pumas_eval.txt"
FILTERED_LIST_DIR="${CODE_DIR_STEP4}/filtered_model_lists"
PUMAS_EVAL_OUTPUT_BASE_DIR="${BASE_PROJECT_DIR}/results/pumas_evaluation_output_filtered"

mkdir -p "${CODE_DIR_STEP4}" "${CODE_DIR_STEP4}/logs" "${WORKER_LOG_DIR}" "${PUMAS_EVAL_OUTPUT_BASE_DIR}"

# --- Check Prerequisites ---
[ -f "${JOB_PARAMS_FILE}" ] || exit 0
[ -d "${FILTERED_LIST_DIR}" ] || exit 0
[ -f "${WORKER_SCRIPT}" ] || exit 0

# --- Submit Array Job ---
actual_job_count=$(wc -l < "${JOB_PARAMS_FILE}")
n_jobs=${1:-$actual_job_count}
[ "$n_jobs" -eq 0 ] && exit 0
[ "$n_jobs" -gt "$actual_job_count" ] && n_jobs=$actual_job_count

sbatch --job-name=PumasEvalFilt \
       --output="${WORKER_LOG_DIR}/PumasEvalFilt_%A_%a.log" \
       --array=1-${n_jobs} \
       "${WORKER_SCRIPT}"
