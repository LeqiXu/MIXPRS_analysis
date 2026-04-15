#!/bin/bash
#SBATCH --job-name=XPASS_submit%j
#SBATCH --output=/home/yd357/pi_paths/pi_zhao/MIXPRS_Revise/Project1_15K_80K/code/4_methods/XPASS/logs/XPASS_submit%j.log
#SBATCH --error=/home/yd357/pi_paths/pi_zhao/MIXPRS_Revise/Project1_15K_80K/code/4_methods/XPASS/logs/XPASS_submit%j.log
#SBATCH --requeue
#SBATCH --mem=1G
#SBATCH --partition=scavenge,day,pi_zhao
#SBATCH -t 0:10:00
#SBATCH -c 1

CODE_DIR="/home/yd357/pi_paths/pi_zhao/MIXPRS_Revise/Project1_15K_80K/code/4_methods/XPASS/codes"
LOG_DIR="/home/yd357/pi_paths/pi_zhao/MIXPRS_Revise/Project1_15K_80K/code/4_methods/XPASS/logs"
job_file="${CODE_DIR}/xpass_params.txt"

mkdir -p "${LOG_DIR}"
> "$job_file"

sim_is=({1..5})
ps=("0.1" "0.01" "0.001" "5e-04")
sample_sizes=("15K" "80K")
pop2s=("EAS" "AFR" "SAS" "AMR")
rhog="0.8"
pop1="EUR"
sample_size_1="80K" # EUR always uses 80K data label

INPUT_DIR="/home/yd357/pi_paths/pi_zhao/MIXPRS_Revise/Project1_15K_80K/data/sim_data/summary_data/discover_validate/XPASS"
RESULT_DIR="/home/yd357/pi_paths/pi_zhao/MIXPRS_Revise/Project1_15K_80K/result/add_benchmark/25.10.4/sim_data/XPASS"

mkdir -p "${RESULT_DIR}"

check_file_exists() {
    local sim_i="$1"
    local p="$2"
    local sample_size_2="$3" 
    local pop2="$4"
    
    local filename="${RESULT_DIR}/sim${sim_i}_p${p}_rho${rhog}_${sample_size_2}_XPASS_EUR_${pop2}_beta_${pop2}.txt"
    if [ -f "$filename" ]; then
        return 0
    else
        return 1
    fi
}

job_count=0
max_jobs=9990 

echo "Scanning for missing XPASS jobs..."

for sim_i in "${sim_is[@]}"; do
    for p in "${ps[@]}"; do
        for sample_size_2 in "${sample_sizes[@]}"; do
            for pop2 in "${pop2s[@]}"; do
                
                if ! check_file_exists "$sim_i" "$p" "$sample_size_2" "$pop2"; then
                    all_inputs_exist=true
                    
                    ss1_file="${INPUT_DIR}/${pop1}_sim${sim_i}_p${p}_rho${rhog}_${sample_size_1}_subcol.txt"
                    if [ ! -f "$ss1_file" ]; then
                        echo "SKIPPING: ss1 file missing: $ss1_file"
                        all_inputs_exist=false
                    fi
                    
                    ss2_file="${INPUT_DIR}/${pop2}_sim${sim_i}_p${p}_rho${rhog}_${sample_size_2}_subcol.txt"
                    if [ ! -f "$ss2_file" ]; then
                        echo "SKIPPING: ss2 file missing: $ss2_file"
                        all_inputs_exist=false
                    fi

                    if "$all_inputs_exist"; then
                        echo "${sim_i},${p},${sample_size_2},${pop2}" >> "$job_file"
                        ((job_count++))
                        
                        if [ $job_count -ge $max_jobs ]; then
                            break 4 
                        fi
                    fi
                fi
            done
            [ $job_count -ge $max_jobs ] && break 3
        done
        [ $job_count -ge $max_jobs ] && break 2
    done
    [ $job_count -ge $max_jobs ] && break
done

job_count=$(wc -l < "$job_file")
if [ "$job_count" -eq 0 ]; then
    echo "No new jobs to submit. All output files already exist."
    exit 0
fi

n_jobs=${1:-$job_count}
sbatch --array=1-$n_jobs ${CODE_DIR}/job_run_xpass.sh
echo "Submitted XPASS array job with $n_jobs tasks."
echo "Parameters written to: $job_file"