#!/bin/bash
#SBATCH --mem=1G
#SBATCH --partition=scavenge,day,week
#SBATCH -t 0:10:00
#SBATCH -c 1

CODE_DIR="/home/yd357/pi_paths/pi_zhao/GWAS_Subsample_JointPRS/v2_PUMAS-EN/PUMAS-EN/codes/step2_subsampling"
GWAS_DIR="/home/yd357/pi_paths/pi_zhao/GWAS_Subsample_JointPRS/v2_PUMAS-EN/PUMAS-EN/data/GWAS"
job_file="${CODE_DIR}/job_params.txt"

> "$job_file"

# Define trait groups and populations
group1=(HDL LDL TC logTG)
group2=(Height BMI SBP DBP PLT)
group3=(WBC NEU LYM MON EOS RBC HCT MCH MCV HB ALT ALP GGT)

# Define population mappings for each group
declare -A group_populations
group_populations[1]="EUR EAS AFR SAS AMR"
group_populations[2]="EUR EAS AFR"
group_populations[3]="EUR EAS"

# Define ensemble methods
ensemble_methods=("EN" "all")

for method in "${ensemble_methods[@]}"; do
    mkdir -p "/home/yd357/pi_paths/pi_zhao/GWAS_Subsample_JointPRS/v2_PUMAS-EN/PUMAS-EN/data/subsampling/PUMA-ensemble.subsampling/${method}"
done

# Generate job parameters by iterating through all combinations
for group_idx in 1 2 3; do
    group_var="group${group_idx}[@]"
    traits=("${!group_var}")
    
    IFS=' ' read -r -a pops <<< "${group_populations[$group_idx]}"
    
    for trait in "${traits[@]}"; do
        for pop in "${pops[@]}"; do
            gwas_file="${trait}_${pop}_inter"

            # Check if any valid GWAS file exists
            if [ -f "${GWAS_DIR}/${gwas_file}.txt" ] || [ -f "${GWAS_DIR}/${gwas_file}.gz" ] || [ -f "${GWAS_DIR}/${gwas_file}.txt.gz" ]; then
                for method in "${ensemble_methods[@]}"; do
                    echo "${trait},${pop},${method}" >> "$job_file"
                done
            else
                echo "Warning: GWAS file not found for ${gwas_file}"
            fi
        done
    done
done

job_count=$(wc -l < "$job_file")
n_jobs=${1:-$job_count}
sbatch --array=1-$n_jobs ${CODE_DIR}/job_run.sh
echo "Submitted PUMA-ensemble subsampling array job with $n_jobs tasks out of $job_count total tasks."
