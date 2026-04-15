#!/bin/bash
#SBATCH --partition=scavenge,day
#SBATCH --requeue
#SBATCH --mem=10G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1 --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --job-name=linear_combine2_array
#SBATCH --output=/gpfs/gibbs/pi/zhao/xz674/logs/IndPRS/1024MIX/linear_combine2_%A_%a.txt

# avoid lmod error
export LMOD_IGNORE_CACHE=1

params_file="/gpfs/gibbs/pi/zhao/xz674/simulation/b.method_calculate/IndPRS/linear_combine/MIXPRS.params.txt"
IFS=' ' read subpop sample1 sample2 h2 sim_i rhog p type type_original model <<< $(sed -n "${SLURM_ARRAY_TASK_ID}p" ${params_file})

if [[ "${sample2}" != "100K" ]]; then
  echo "Please provide the available sample size as 100K."
fi

if [[ "${sample1}" == "UKB" ]]; then
  sample1="ukbb"
  weight_sample1="UKB"
else
  weight_sample1=${sample1}
fi

JointPRS_EUR="/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/JointPRS/sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${sample1}_${sample2}_JointPRS_real_EUR_EAS_AFR_SAS_AMR_beta_EUR.txt"
JointPRS_EAS="/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/JointPRS/sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${sample1}_${sample2}_JointPRS_real_EUR_EAS_AFR_SAS_AMR_beta_EAS.txt"
JointPRS_AFR="/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/JointPRS/sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${sample1}_${sample2}_JointPRS_real_EUR_EAS_AFR_SAS_AMR_beta_AFR.txt"
JointPRS_SAS="/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/JointPRS/sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${sample1}_${sample2}_JointPRS_real_EUR_EAS_AFR_SAS_AMR_beta_SAS.txt"
JointPRS_AMR="/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/JointPRS/sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${sample1}_${sample2}_JointPRS_real_EUR_EAS_AFR_SAS_AMR_beta_AMR.txt"

if [[ "${subpop}" == "EUR" ]]; then
    SDPRX_EUR="/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/SDPRX/sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${sample1}_${sample2}_SDPRX_real_EUR_EAS_beta_EUR.txt"
    sample_size="${sample1}"
else
    SDPRX_EUR="/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/SDPRX/sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${sample1}_${sample2}_SDPRX_real_EUR_${subpop}_beta_EUR.txt"
    sample_size="${sample2}"
fi

SDPRX_EAS="/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/SDPRX/sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${sample1}_${sample2}_SDPRX_real_EUR_EAS_beta_EAS.txt"
SDPRX_AFR="/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/SDPRX/sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${sample1}_${sample2}_SDPRX_real_EUR_AFR_beta_AFR.txt"
SDPRX_SAS="/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/SDPRX/sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${sample1}_${sample2}_SDPRX_real_EUR_SAS_beta_SAS.txt"
SDPRX_AMR="/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/SDPRX/sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${sample1}_${sample2}_SDPRX_real_EUR_AMR_beta_AMR.txt"

prs_file="${JointPRS_EUR},${JointPRS_EAS},${JointPRS_AFR},${JointPRS_SAS},${JointPRS_AMR},${SDPRX_EUR},${SDPRX_EAS},${SDPRX_AFR},${SDPRX_SAS},${SDPRX_AMR}"

weight_file1="/gpfs/gibbs/pi/zhao/xz674/result/sim_result/IndPRS/Final_weight/sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_EUR_EAS_AFR_SAS_AMR_${weight_sample1}_${sample2}_JointPRS_SDPRX_${type}_fold1_${subpop}_${model}_weights.txt"
weight_file2="/gpfs/gibbs/pi/zhao/xz674/result/sim_result/IndPRS/Final_weight/sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_EUR_EAS_AFR_SAS_AMR_${weight_sample1}_${sample2}_JointPRS_SDPRX_${type}_fold2_${subpop}_${model}_weights.txt"
weight_file3="/gpfs/gibbs/pi/zhao/xz674/result/sim_result/IndPRS/Final_weight/sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_EUR_EAS_AFR_SAS_AMR_${weight_sample1}_${sample2}_JointPRS_SDPRX_${type}_fold3_${subpop}_${model}_weights.txt"
weight_file4="/gpfs/gibbs/pi/zhao/xz674/result/sim_result/IndPRS/Final_weight/sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_EUR_EAS_AFR_SAS_AMR_${weight_sample1}_${sample2}_JointPRS_SDPRX_${type}_fold4_${subpop}_${model}_weights.txt"
weight_file="${weight_file1},${weight_file2},${weight_file3},${weight_file4}"

sst_file="/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/sim_data/summary_data/discover_validate/MIX/${subpop}_sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${sample_size}_MIX_real_all.txt"

out_dir="/gpfs/gibbs/pi/zhao/xz674/result/sim_result/IndPRS/MIXPRS"
mkdir -p ${out_dir}
out_name="sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_EUR_EAS_AFR_SAS_AMR_${sample1}_${sample2}_JointPRS_SDPRX_${type}_${model}_weights"

output_file="${out_dir}/${out_name}_${subpop}_MIXPRS.txt"

if [[ ! -e ${output_file} ]]; then

module load miniconda
conda activate py_env

if [[ "${subpop}" == "EUR" ]]; then

python /gpfs/gibbs/pi/zhao/xz674/MIXPRS/main/MIX_final_combine.py \
--ref_dir=/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ref_data/PRScsx/ukbb \
--sst_file=${sst_file} \
--pop=${subpop} \
--prs_beta_file=${prs_file} \
--weight_file=${weight_file} \
--indep_approx=${approx} \
--out_dir=${out_dir} \
--out_name=${out_name}

else

python /gpfs/gibbs/pi/zhao/xz674/MIXPRS/main/MIX_final_combine.py \
--ref_dir=/gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/ref_data/discover_validate/100K \
--sst_file=${sst_file} \
--pop=${subpop} \
--prs_beta_file=${prs_file} \
--weight_file=${weight_file} \
--indep_approx=${approx} \
--out_dir=${out_dir} \
--out_name=${out_name}

fi

fi

