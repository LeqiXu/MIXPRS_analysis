#!/bin/bash
#SBATCH --job-name=JointPRS_submit%j
#SBATCH --output=/home/yd357/pi_paths/pi_zhao/MIXPRS_Revise/Project1_15K_80K/code/4_methods/JointPRS/logs/25.10.5/JointPRS_submit%j.log
#SBATCH --error=/home/yd357/pi_paths/pi_zhao/MIXPRS_Revise/Project1_15K_80K/code/4_methods/JointPRS/logs/25.10.5/JointPRS_submit%j.log
#SBATCH --requeue
#SBATCH --mem=1G
#SBATCH --partition=scavenge,day,pi_zhao
#SBATCH -t 0:10:00
#SBATCH -c 1

CODE_DIR="/home/yd357/pi_paths/pi_zhao/MIXPRS_Revise/Project1_15K_80K/code/4_methods/JointPRS/codes"
job_file="${CODE_DIR}/jointprs_params.txt"

> "$job_file"

pops=("EUR" "EAS" "AFR" "SAS" "AMR")
sim_is=({1..5})
ps=("0.1" "0.01" "0.001" "5e-04")
sample_sizes=("15K" "80K")
rhog="0.8"
chrs=({1..22})

INPUT_DIR="/home/yd357/pi_paths/pi_zhao/MIXPRS_Revise/Project1_15K_80K/data/sim_data/summary_data/discover_validate/JointPRS"
REF_DIR="/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ref_data/PRScsx/1kg"
BIM_PREFIX="/gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/data/sim_data/geno_data/All/All_test"
RESULT_DIR="/home/yd357/pi_paths/pi_zhao/MIXPRS_Revise/Project1_15K_80K/result/add_benchmark/25.10.4/sim_data/JointPRS"

mkdir -p "${RESULT_DIR}"

check_file_exists() {
    local sim_i="$1"
    local p="$2"
    local sample_size="$3"
    local chr="$4"
    local pop="$5"
    
    local filename="${RESULT_DIR}/sim${sim_i}_p${p}_rho${rhog}_${sample_size}_JointPRS_${pop}_pst_eff_a1_b0.5_phiauto_chr${chr}.txt"
    if [ -f "$filename" ]; then
        return 0 
    else
        return 1  
    fi
}

job_count=0
max_jobs=9990

for sim_i in "${sim_is[@]}"; do
    for p in "${ps[@]}"; do
        for sample_size in "${sample_sizes[@]}"; do
            for chr in "${chrs[@]}"; do
                if ! check_file_exists "$sim_i" "$p" "$sample_size" "$chr" "EAS"; then
                    
                    all_files_exist=true
                    
                    # 1. Check for the EUR file (always using 80K)
                    eur_input_file="${INPUT_DIR}/EUR_sim${sim_i}_p${p}_rho${rhog}_80K_subcol.txt"
                    if [ ! -f "$eur_input_file" ]; then
                        all_files_exist=false
                    fi
                    
                    # 2. Check for other population files (using the current loop sample_size)
                    if "$all_files_exist"; then
                        for pop in "EAS" "AFR" "SAS" "AMR"; do
                            input_file="${INPUT_DIR}/${pop}_sim${sim_i}_p${p}_rho${rhog}_${sample_size}_subcol.txt"
                            if [ ! -f "$input_file" ]; then
                                all_files_exist=false
                                break
                            fi
                        done
                    fi
                    
                    if "$all_files_exist"; then
                        echo "${sim_i},${p},${sample_size},${chr}" >> "$job_file"
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
n_jobs=${1:-$job_count}
sbatch --array=1-$n_jobs ${CODE_DIR}/job_run_jointprs.sh
echo "Submitted JointPRS array job with $n_jobs tasks out of $job_count total tasks."
echo "Parameters written to: $job_file"