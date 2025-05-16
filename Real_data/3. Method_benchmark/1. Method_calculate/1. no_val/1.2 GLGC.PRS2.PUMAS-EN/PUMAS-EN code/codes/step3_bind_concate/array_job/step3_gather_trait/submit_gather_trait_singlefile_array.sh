#!/bin/bash
#SBATCH --mem=1G
#SBATCH --partition=scavenge,day,week,pi_zhao
#SBATCH -t 0:20:00
#SBATCH -c 1

# --- Configuration ---
BASE_PROJECT_DIR="/home/yd357/pi_paths/pi_zhao/GWAS_Subsample_JointPRS/v3/concatenate"
CODE_DIR_STEP3="${BASE_PROJECT_DIR}/codes/array_job/step3_gather_trait"
SUBMITTER_LOG_DIR="${BASE_PROJECT_DIR}/logs/step3_submitter"
WORKER_LOG_DIR="${BASE_PROJECT_DIR}/logs/step3_gather_trait_worker"
JOB_PARAMS_FILE="${CODE_DIR_STEP3}/job_params_step3_singlefile.txt"
POP_INFO_FILE="/home/yd357/pi_paths/pi_zhao/GWAS_Subsample_JointPRS/v2_PUMAS-EN/PUMAS-EN/data/basic_pop_info/basic_pop_info.csv"
FINAL_OUTPUT_BASE_DIR="${BASE_PROJECT_DIR}/results/concatenation_trait_single_file"

mkdir -p "${CODE_DIR_STEP3}" "${SUBMITTER_LOG_DIR}" "${WORKER_LOG_DIR}" "${FINAL_OUTPUT_BASE_DIR}"
> "${JOB_PARAMS_FILE}"

# --- Generate Job Parameters ---
CASES=("PUMASite1" "PUMASite2" "PUMASite3" "PUMASite4" "Regular")

tail -n +2 "${POP_INFO_FILE}" | while IFS=, read -r trait ancestry _
do
    target_trait=$(echo "$trait" | tr -d '"[:space:]')
    target_ancestry=$(echo "$ancestry" | tr -d '"[:space:]')

    [[ -z "$target_trait" || -z "$target_ancestry" ]] && continue

    for case_type in "${CASES[@]}"; do
        echo "${target_trait},${target_ancestry},${case_type}" >> "${JOB_PARAMS_FILE}"
    done
done

# --- Submit SLURM Array Job ---
actual_job_count=$(wc -l < "${JOB_PARAMS_FILE}")
n_jobs=${1:-$actual_job_count}
n_jobs=$(( n_jobs > actual_job_count ? actual_job_count : n_jobs ))

if [ "$n_jobs" -gt 0 ]; then
    sbatch --job-name=GatherSingle \
           --output="${WORKER_LOG_DIR}/GatherSingle_%A_%a.log" \
           --array=1-$n_jobs \
           "${CODE_DIR_STEP3}/job_run_gather_trait_singlefile.sh"
fi
