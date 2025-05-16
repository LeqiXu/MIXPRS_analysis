#!/bin/bash
#SBATCH --mem=1G
#SBATCH --partition=scavenge,day,week
#SBATCH -t 0:10:00
#SBATCH -c 1

CODE_DIR="/home/yd357/pi_paths/pi_zhao/GWAS_Subsample_JointPRS/v2_PUMAS-EN/PUMAS-EN/codes/step1_gwas_qc/1.frq_transforming"
job_file="${CODE_DIR}/job_params.txt"

> "$job_file"

# Define populations to process
pops=("EUR" "AMR" "EAS" "SAS" "AFR")

# Add each population to the parameter file
for pop in "${pops[@]}"; do
    echo "${pop}" >> "$job_file"
done

# Count total jobs and submit array
job_count=$(wc -l < "$job_file")
n_jobs=${1:-$job_count}
sbatch --array=1-$n_jobs ${CODE_DIR}/job_run.sh
echo "Submitted frequency conversion array job with $n_jobs tasks out of $job_count total tasks."