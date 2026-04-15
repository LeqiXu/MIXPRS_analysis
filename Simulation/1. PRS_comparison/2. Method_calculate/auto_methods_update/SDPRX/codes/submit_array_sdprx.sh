#!/bin/bash
#SBATCH --job-name=SDPRX_submit%j
#SBATCH --output=/home/yd357/pi_paths/pi_zhao/MIXPRS_Revise/Project1_15K_80K/code/4_methods/SDPRX/logs/11.2/SDPRX_submit%j.log
#SBATCH --error=/home/yd357/pi_paths/pi_zhao/MIXPRS_Revise/Project1_15K_80K/code/4_methods/SDPRX/logs/11.2/SDPRX_submit%j.log
#SBATCH --requeue
#SBATCH --mem=1G
#SBATCH --partition=scavenge,day,pi_zhao
#SBATCH -t 0:10:00
#SBATCH -c 1

CODE_DIR="/home/yd357/pi_paths/pi_zhao/MIXPRS_Revise/Project1_15K_80K/code/4_methods/SDPRX/codes"
job_file="${CODE_DIR}/sdprx_params.txt"

> "$job_file"

# --- Parameters ---
sim_is=({1..5})
ps=("0.1" "0.01" "0.001" "5e-04")
sample_sizes=("15K" "80K") # This is for sample 2 (pop2)
pop2s=("EAS" "AFR" "SAS" "AMR") # Paired with EUR
rhog="0.8"
chrs=({1..22})
pop1="EUR"
sample_size_1="80K" # EUR always uses 80K data label

# --- Paths ---
# Assumes SDPRX-formatted summary stats are in a parallel dir to JointPRS
INPUT_DIR="/home/yd357/pi_paths/pi_zhao/MIXPRS_Revise/Project1_15K_80K/data/sim_data/summary_data/discover_validate/SDPRX"
RESULT_DIR="/home/yd357/pi_paths/pi_zhao/MIXPRS_Revise/Project1_15K_80K/result/add_benchmark/25.10.4/sim_data/SDPRX"
# Assumes popcorn results are in a parallel dir
POPCORN_DIR="/home/yd357/pi_paths/pi_zhao/MIXPRS_Revise/Project1_15K_80K/result/add_benchmark/25.10.4/sim_data/popcorn"

mkdir -p "${RESULT_DIR}"
mkdir -p "/home/yd357/pi_paths/pi_zhao/MIXPRS_Revise/Project1_15K_80K/code/4_methods/SDPRX/logs"

# Function to check if the final output file exists
check_file_exists() {
    local sim_i="$1"
    local p="$2"
    local sample_size_2="$3" # 15K or 80K
    local pop2="$4"
    local chr="$5"
    
    # SDPRX appends _1.txt and _2.txt. We check for _2.txt as the sentinel.
    local filename="${RESULT_DIR}/sim${sim_i}_p${p}_rho${rhog}_${sample_size_2}_${pop1}_${pop2}_SDPRX_chr${chr}_2.txt"
    if [ -f "$filename" ]; then
        return 0 # File exists
    else
        return 1 # File does not exist
    fi
}

job_count=0
max_jobs=9990

echo "Scanning for missing jobs..."

for sim_i in "${sim_is[@]}"; do
    for p in "${ps[@]}"; do
        for sample_size_2 in "${sample_sizes[@]}"; do
            for pop2 in "${pop2s[@]}"; do
                
                # Check for required input files for this (sim_i, p, sample_size_2, pop2) combo
                all_inputs_exist=true
                
                # 1. Check Popcorn file (provides rho_est)
                popcorn_file="${POPCORN_DIR}/sim${sim_i}_p${p}_rho${rhog}_${pop1}_${pop2}_${sample_size_1}_${sample_size_2}_popcorn_corr.txt"
                if [ ! -f "$popcorn_file" ]; then
                    echo "SKIPPING: Popcorn file missing for $popcorn_file"
                    all_inputs_exist=false
                fi
                
                # 2. Check ss1 (EUR) file
                ss1_file="${INPUT_DIR}/${pop1}_sim${sim_i}_p${p}_rho${rhog}_${sample_size_1}_subcol.txt"
                if [ ! -f "$ss1_file" ]; then
                    echo "SKIPPING: ss1 file missing for $ss1_file"
                    all_inputs_exist=false
                fi
                
                # 3. Check ss2 (pop2) file
                ss2_file="${INPUT_DIR}/${pop2}_sim${sim_i}_p${p}_rho${rhog}_${sample_size_2}_subcol.txt"
                if [ ! -f "$ss2_file" ]; then
                    echo "SKIPPING: ss2 file missing for $ss2_file"
                    all_inputs_exist=false
                fi

                if "$all_inputs_exist"; then
                    for chr in "${chrs[@]}"; do
                        if ! check_file_exists "$sim_i" "$p" "$sample_size_2" "$pop2" "$chr"; then
                            
                            echo "${sim_i},${p},${sample_size_2},${pop2},${chr}" >> "$job_file"
                            ((job_count++))
                            
                            if [ $job_count -ge $max_jobs ]; then
                                break 5 # Break all loops
                            fi
                        fi
                    done
                fi
            done
            [ $job_count -ge $max_jobs ] && break 3
        done
        [ $job_count -ge $max_jobs ] && break 2
    done
    [ $job_count -ge $max_jobs ] && break
done

job_count=$(wc -l < "$job_file")
n_jobs=${1:-$job_count}
sbatch --array=1-$n_jobs ${CODE_DIR}/job_run_sdprx.sh
echo "Submitted SDPRX array job with $n_jobs tasks."
echo "Parameters written to: $job_file"