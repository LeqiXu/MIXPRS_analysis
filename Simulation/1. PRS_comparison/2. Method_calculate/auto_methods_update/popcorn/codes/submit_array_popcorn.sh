#!/bin/bash
#SBATCH --job-name=Popcorn_submit%j
#SBATCH --output=/home/yd357/pi_paths/pi_zhao/MIXPRS_Revise/Project1_15K_80K/code/4_methods/popcorn/logs/25.10.29/Popcorn_submit%j.log
#SBATCH --error=/home/yd357/pi_paths/pi_zhao/MIXPRS_Revise/Project1_15K_80K/code/4_methods/popcorn/logs/25.10.29/Popcorn_submit%j.log
#SBATCH --requeue
#SBATCH --mem=1G
#SBATCH --partition=scavenge,day,pi_zhao
#SBATCH -t 0:10:00
#SBATCH -c 1

CODE_DIR="/home/yd357/pi_paths/pi_zhao/MIXPRS_Revise/Project1_15K_80K/code/4_methods/popcorn/codes"
LOG_DIR="/home/yd357/pi_paths/pi_zhao/MIXPRS_Revise/Project1_15K_80K/code/4_methods/popcorn/logs/25.10.29"
job_file="${CODE_DIR}/popcorn_params.txt"

> "$job_file"
 
# --- Parameters ---
sim_is=({1..5})
ps=("0.1" "0.01" "0.001" "5e-04")
sample_sizes=("15K" "80K") # Sample sizes for non-EUR pops
pop2s=("EAS" "AFR" "SAS" "AMR") # Target populations (pop1 is EUR)
rhog="0.8"
pop1="EUR"
sample1_label="80K" # EUR always uses 80K data

# --- Directories ---
INPUT_DIR="/home/yd357/pi_paths/pi_zhao/MIXPRS_Revise/Project1_15K_80K/data/sim_data/summary_data/discover_validate/popcorn"
RESULT_DIR="/home/yd357/pi_paths/pi_zhao/MIXPRS_Revise/Project1_15K_80K/result/add_benchmark/25.10.4/sim_data/popcorn"

mkdir -p "${RESULT_DIR}"
mkdir -p "${LOG_DIR}"

job_count=0
max_jobs=9990

echo "Generating parameter file at $job_file..."

for sim_i in "${sim_is[@]}"; do
    for p in "${ps[@]}"; do
        for sample_size in "${sample_sizes[@]}"; do 
            for pop2 in "${pop2s[@]}"; do
                
                output_file="${RESULT_DIR}/sim${sim_i}_p${p}_rho${rhog}_${pop1}_${pop2}_${sample1_label}_${sample_size}_popcorn_corr.txt"
                input_file1="${INPUT_DIR}/${pop1}_sim${sim_i}_p${p}_rho${rhog}_${sample1_label}_subcol.txt"
                input_file2="${INPUT_DIR}/${pop2}_sim${sim_i}_p${p}_rho${rhog}_${sample_size}_subcol.txt"

                if [ -f "$output_file" ]; then
                    continue
                fi

                if [ ! -f "$input_file1" ] || [ ! -f "$input_file2" ]; then
                    echo "Warning: Missing input files for sim $sim_i, p $p, size $sample_size, pop $pop2. Skipping."
                    [ ! -f "$input_file1" ] && echo "  Missing: $input_file1"
                    [ ! -f "$input_file2" ] && echo "  Missing: $input_file2"
                    continue
                fi

                echo "${sim_i},${p},${sample_size},${pop2}" >> "$job_file"
                ((job_count++))
                
                if [ $job_count -ge $max_jobs ]; then
                    break 4
                fi
            done
            [ $job_count -ge $max_jobs ] && break 3
        done
        [ $job_count -ge $max_jobs ] && break 2
    done
    [ $job_count -ge $max_jobs ] && break
done

n_jobs=${1:-$job_count}

if [ "$n_jobs" -eq 0 ]; then
    echo "No new jobs to submit. All output files may already exist or input files are missing."
else
    sbatch --array=1-$n_jobs ${CODE_DIR}/job_run_popcorn.sh
    echo "Submitted Popcorn array job with $n_jobs tasks (out of $job_count total runnable tasks)."
    echo "Parameters written to: $job_file"
fi