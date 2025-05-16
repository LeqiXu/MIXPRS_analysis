#!/bin/bash
#SBATCH --mem=1G
#SBATCH --partition=scavenge,pi_zhao
#SBATCH -t 0:05:00
#SBATCH -c 1

BASE_DIR="/home/yd357/pi_paths/pi_zhao/GWAS_Subsample_JointPRS/v2_PUMAS-EN/single-pop/PRS-CS"
CODE_DIR="${BASE_DIR}/codes/base" 
LOG_DIR="${BASE_DIR}/logs/formatting"
PARAM_DIR="${CODE_DIR}"
job_param_file="${PARAM_DIR}/format_params.txt"
job_script="${CODE_DIR}/format_gwas_job.sh"

mkdir -p "${LOG_DIR}"
mkdir -p "${PARAM_DIR}"

> "$job_param_file"
echo "Cleared previous parameter file: ${job_param_file}"

# Read population and trait information
POP_INFO="/home/yd357/pi_paths/pi_zhao/GWAS_Subsample_JointPRS/v2_PUMAS-EN/PUMAS-EN/data/basic_pop_info/basic_pop_info.csv"
if [ ! -f "$POP_INFO" ]; then
    echo "ERROR: Population info file not found: ${POP_INFO}"
    exit 1
fi

echo "Reading population info from: ${POP_INFO}"
processed_count=0
skipped_count=0

tail -n +2 "${POP_INFO}" | while IFS=, read -r trait ancestry N h2; do
    trait=$(echo $trait | tr -d '"')
    ancestry=$(echo $ancestry | tr -d '"')

    # --- Process PUMAS=T mode (with iterations) ---
    pumas="T"
    gwas_base_path="/home/yd357/pi_paths/pi_zhao/GWAS_Subsample_JointPRS/v2_PUMAS-EN/PUMAS-EN/data/subsampling/PUMA-ensemble.subsampling/EN/${trait}_${ancestry}_inter"

    if [ -d "$gwas_base_path" ]; then
        for iter in $(seq 1 4); do
            gwas_raw_file="${gwas_base_path}/${trait}_${ancestry}_inter.gwas.omnibus.ite${iter}.txt"

            if [ -f "$gwas_raw_file" ]; then
                 echo "${trait},${ancestry},${pumas},${iter}" >> "$job_param_file"
                 ((processed_count++))
            else
                 echo "INFO: Skipping PUMAS=T - Raw GWAS file not found: ${gwas_raw_file}"
                 ((skipped_count++))
            fi
        done
    else
        echo "INFO: Skipping PUMAS=T - Base GWAS directory not found for ${trait}_${ancestry}: ${gwas_base_path}"
        ((skipped_count+=4)) # Skip all 4 iterations if base dir missing
    fi

    # --- Process PUMAS=F mode (no iterations) ---
    pumas="F"
    gwas_raw_file="/gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/data/summary_data/clean/${trait}_${ancestry}_inter_clean.txt"

    if [ -f "$gwas_raw_file" ]; then
         echo "${trait},${ancestry},${pumas}," >> "$job_param_file"
         ((processed_count++))
    else
         echo "INFO: Skipping PUMAS=F - Raw GWAS file not found: ${gwas_raw_file}"
         ((skipped_count++))
    fi

done

echo "Finished generating parameter list."
echo "Found $processed_count potential formatting tasks."
echo "Skipped $skipped_count tasks due to missing raw files or directories."

job_count=$(wc -l < "$job_param_file")
max_jobs=9990

if [ $job_count -gt 0 ]; then
    echo "Submitting array job with $job_count tasks..."
    sbatch --array=1-$job_count ${job_script}
    echo "Submitted formatting array job with $job_count tasks."
else
    echo "No formatting jobs to submit. All raw input files might be missing or no combinations found."
fi

echo "Submission script finished."