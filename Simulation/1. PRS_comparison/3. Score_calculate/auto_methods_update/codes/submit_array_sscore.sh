#!/bin/bash
#SBATCH --job-name=sscore_submit%j
#SBATCH --output=/home/yd357/pi_paths/pi_zhao/MIXPRS_Revise/Project1_15K_80K/code/5_evaluation/1_sscore/logs/25.11.8.1/sscore_submit%j.log
#SBATCH --requeue
#SBATCH --mem=1G
#SBATCH --partition=scavenge,day,pi_zhao
#SBATCH -t 0:10:00
#SBATCH -c 1

CODE_DIR="/home/yd357/pi_paths/pi_zhao/MIXPRS_Revise/Project1_15K_80K/code/5_evaluation/1_sscore/codes"
job_file="${CODE_DIR}/sscore_params.txt"

> "$job_file"

# --- Parameters ---
sim_is=({1..5})
ps=("0.1" "0.01" "0.001" "5e-04")
sample_sizes=("15K" "80K") 
pop2s=("EAS" "AFR" "SAS" "AMR")
rhog="0.8"
pop1="EUR"

# --- Score Output Directories ---
SCORE_DIR_BASE="/home/yd357/pi_paths/pi_zhao/MIXPRS_Revise/Project1_15K_80K/result/add_benchmark/25.10.4/sim_data/prs_scores"
SCORE_DIR_JOINTPRS="${SCORE_DIR_BASE}/JointPRS"
SCORE_DIR_SDPRX="${SCORE_DIR_BASE}/SDPRX"
SCORE_DIR_XPASS="${SCORE_DIR_BASE}/XPASS"
SCORE_DIR_MIXPRS_J="${SCORE_DIR_BASE}/MIXPRS/MIXPRS-JointPRS"
SCORE_DIR_MIXPRS_S="${SCORE_DIR_BASE}/MIXPRS/MIXPRS-SDPRX"

mkdir -p "${SCORE_DIR_JOINTPRS}" "${SCORE_DIR_SDPRX}" "${SCORE_DIR_XPASS}" \
           "${SCORE_DIR_MIXPRS_J}" "${SCORE_DIR_MIXPRS_S}"

job_count=0
echo "Scanning for missing .sscore files..."

for sim_i in "${sim_is[@]}"; do
  for p in "${ps[@]}"; do
    for sample_size in "${sample_sizes[@]}"; do
      for pop2 in "${pop2s[@]}"; do
        
        # --- 1. Check JointPRS (Original) ---
        sentinel_file="${SCORE_DIR_JOINTPRS}/sim${sim_i}_p${p}_rho${rhog}_${sample_size}_JointPRS_beta_AMR_prs_${pop2}.sscore"
        if [ ! -f "$sentinel_file" ]; then
          echo "JointPRS,${sim_i},${p},${sample_size},${pop2}" >> "$job_file"
          ((job_count++))
        fi

        # --- 2. Check SDPRX (Original) ---
        sentinel_file="${SCORE_DIR_SDPRX}/sim${sim_i}_p${p}_rho${rhog}_${sample_size}_${pop1}_${pop2}_beta_${pop2}_prs_${pop2}.sscore"
        if [ ! -f "$sentinel_file" ]; then
          echo "SDPRX,${sim_i},${p},${sample_size},${pop2}" >> "$job_file"
          ((job_count++))
        fi
        
        # --- 3. Check XPASS ---
        sentinel_file="${SCORE_DIR_XPASS}/sim${sim_i}_p${p}_rho${rhog}_${sample_size}_XPASS_${pop1}_${pop2}_beta_${pop2}_prs_${pop2}.sscore"
        if [ ! -f "$sentinel_file" ]; then
          echo "XPASS,${sim_i},${p},${sample_size},${pop2}" >> "$job_file"
          ((job_count++))
        fi


      done
    done
  done
done

job_count=$(wc -l < "$job_file")
if [ "$job_count" -eq 0 ]; then
    echo "No new .sscore jobs to submit. All files already exist."
    exit 0
fi

n_jobs=${1:-$job_count}
# Submit with a throttle (e.g., 100 jobs at a time)
sbatch --array=1-$n_jobs%100 ${CODE_DIR}/job_run_sscore.sh
echo "Submitted .sscore array job with $n_jobs tasks."
echo "Parameters written to: $job_file"