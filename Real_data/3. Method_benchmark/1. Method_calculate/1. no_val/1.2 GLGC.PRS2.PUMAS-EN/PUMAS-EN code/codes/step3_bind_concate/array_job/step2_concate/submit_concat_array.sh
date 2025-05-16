#!/bin/bash
#SBATCH --mem=1G
#SBATCH --partition=scavenge,day,week,pi_zhao
#SBATCH -t 0:30:00
#SBATCH -c 1

# --- Configuration ---
BASE_DIR="/home/yd357/pi_paths/pi_zhao/GWAS_Subsample_JointPRS/v3/concatenate"
CODE_DIR="${BASE_DIR}/codes/array_job"
LOG_DIR="${BASE_DIR}/codes/logs"
JOB_PARAMS_FILE="${CODE_DIR}/job_params.txt"
POP_INFO="/home/yd357/pi_paths/pi_zhao/GWAS_Subsample_JointPRS/v2_PUMAS-EN/PUMAS-EN/data/basic_pop_info/basic_pop_info.csv"

CONCAT_SINGLE_POP_DIR="${BASE_DIR}/results"
CONCAT_FINAL_DIR="${BASE_DIR}/results/concatenation"

mkdir -p "${CONCAT_SINGLE_POP_DIR}" "${CONCAT_FINAL_DIR}" "${LOG_DIR}"
> "${JOB_PARAMS_FILE}"

# --- Generate Job Parameters ---
CASES=("PUMASite1" "PUMASite2" "PUMASite3" "PUMASite4" "Regular")
tail -n +2 "${POP_INFO}" | while IFS=, read -r trait ancestry N h2 || [[ -n "$trait" ]]; do
    trait=$(echo "$trait" | tr -d '"[:space:]')
    ancestry=$(echo "$ancestry" | tr -d '"[:space:]')
    [[ -z "$trait" || -z "$ancestry" ]] && continue
    for case in "${CASES[@]}"; do
        echo "${trait},${ancestry},${case}" >> "${JOB_PARAMS_FILE}"
    done
done

job_count=$(wc -l < "${JOB_PARAMS_FILE}")
n_jobs=${1:-$job_count}
[[ $n_jobs -gt $job_count ]] && n_jobs=$job_count

# --- Submit Job ---
if [ $n_jobs -gt 0 ]; then
    sbatch --job-name=Concat4method \
           --output="${LOG_DIR}/Concat4method_%A_%a.log" \
           --array=1-$n_jobs \
           "${CODE_DIR}/job_run_concat.sh"
fi
