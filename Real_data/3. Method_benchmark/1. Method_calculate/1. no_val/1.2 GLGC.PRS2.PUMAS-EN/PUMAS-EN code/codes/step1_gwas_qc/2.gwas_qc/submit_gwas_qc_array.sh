#!/bin/bash
#SBATCH --mem=1G
#SBATCH --partition=scavenge,day,week
#SBATCH -t 0:10:00
#SBATCH -c 1

CODE_DIR="/home/yd357/pi_paths/pi_zhao/GWAS_Subsample_JointPRS/v2_PUMAS-EN/PUMAS-EN/codes/step1_gwas_qc/2.gwas_qc"
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

declare -A group_names
group_names[1]="GLGC"
group_names[2]="PAGE"
group_names[3]="BBJ"

# Sample sizes
declare -A sample_sizes
# Group 1
sample_sizes[HDL_EUR]=885546; sample_sizes[HDL_EAS]=116404; sample_sizes[HDL_AFR]=90804; sample_sizes[HDL_SAS]=33953; sample_sizes[HDL_AMR]=47276
sample_sizes[LDL_EUR]=840012; sample_sizes[LDL_EAS]=79693; sample_sizes[LDL_AFR]=87759; sample_sizes[LDL_SAS]=33658; sample_sizes[LDL_AMR]=33989
sample_sizes[TC_EUR]=929739; sample_sizes[TC_EAS]=144579; sample_sizes[TC_AFR]=92554; sample_sizes[TC_SAS]=34135; sample_sizes[TC_AMR]=48055
sample_sizes[logTG_EUR]=860679; sample_sizes[logTG_EAS]=81071; sample_sizes[logTG_AFR]=89467; sample_sizes[logTG_SAS]=34023; sample_sizes[logTG_AMR]=37273
# Group 2
sample_sizes[Height_EUR]=252357; sample_sizes[Height_EAS]=159095; sample_sizes[Height_AFR]=49781
sample_sizes[BMI_EUR]=233787; sample_sizes[BMI_EAS]=158284; sample_sizes[BMI_AFR]=49335
sample_sizes[SBP_EUR]=728893; sample_sizes[SBP_EAS]=179000; sample_sizes[SBP_AFR]=35433
sample_sizes[DBP_EUR]=746038; sample_sizes[DBP_EAS]=179000; sample_sizes[DBP_AFR]=35433
sample_sizes[PLT_EUR]=539667; sample_sizes[PLT_EAS]=179000; sample_sizes[PLT_AFR]=29328
# Group 3
sample_sizes[WBC_EUR]=559083; sample_sizes[WBC_EAS]=179000
sample_sizes[NEU_EUR]=517889; sample_sizes[NEU_EAS]=179000
sample_sizes[LYM_EUR]=523524; sample_sizes[LYM_EAS]=179000
sample_sizes[MON_EUR]=520195; sample_sizes[MON_EAS]=179000
sample_sizes[EOS_EUR]=473152; sample_sizes[EOS_EAS]=179000
sample_sizes[RBC_EUR]=542043; sample_sizes[RBC_EAS]=179000
sample_sizes[HCT_EUR]=559099; sample_sizes[HCT_EAS]=179000
sample_sizes[MCH_EUR]=483664; sample_sizes[MCH_EAS]=179000
sample_sizes[MCV_EUR]=540967; sample_sizes[MCV_EAS]=179000
sample_sizes[HB_EUR]=408112; sample_sizes[HB_EAS]=179000
sample_sizes[ALT_EUR]=437267; sample_sizes[ALT_EAS]=179000
sample_sizes[ALP_EUR]=437267; sample_sizes[ALP_EAS]=179000
sample_sizes[GGT_EUR]=437267; sample_sizes[GGT_EAS]=179000

# Generate job parameters by iterating through all combinations
for group_idx in 1 2 3; do
    group_var="group${group_idx}[@]"
    traits=("${!group_var}")
    
    IFS=' ' read -r -a pops <<< "${group_populations[$group_idx]}"
    
    for trait in "${traits[@]}"; do
        for pop in "${pops[@]}"; do
            # Get the sample size for this trait-population combination
            n_samples=${sample_sizes[${trait}_${pop}]}
            group_name=${group_names[$group_idx]}

            # Check if input file exists
            input_file="/home/yd357/pi_paths/pi_zhao/GWAS_Subsample_JointPRS/data/summary_data/X_Wing/${trait}_${pop}_inter.txt"


            if [ -f "$input_file" ]; then
                echo "${trait},${pop},${n_samples},${group_name}" >> "$job_file"
            else
                echo "Warning: Input file not found: ${input_file}"
            fi
        done
    done
done

# Count total jobs and submit array
job_count=$(wc -l < "$job_file")
n_jobs=${1:-$job_count}
sbatch --array=1-$n_jobs ${CODE_DIR}/job_run.sh
echo "Submitted GWAS QC array job with $n_jobs tasks out of $job_count total tasks."