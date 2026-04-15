#!/bin/bash
#SBATCH --partition=scavenge,day
#SBATCH --requeue
#SBATCH --mem=15G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1 --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --job-name=SDPRX_calc_array
#SBATCH --output=/gpfs/gibbs/pi/zhao/xz674/logs/IndPRS/1013SDPRX_calc/SDPRX_calc_%A_%a.txt

# avoid lmod error
export LMOD_IGNORE_CACHE=1

params_file="/gpfs/gibbs/pi/zhao/xz674/simulation/b.method_calculate/IndPRS/SDPRX/params.txt"
IFS=' ' read pop1 sample1 pop2 sample2 subpop h2 sim_i rhog p ff chr <<< $(sed -n "${SLURM_ARRAY_TASK_ID}p" ${params_file})

if [[ ${sample1} == "UKB" ]]; then
if [[ ${sample2} == "100K" ]]; then
sample1_popcorn="ukbb"
fi
else
sample1_popcorn=${sample1}
fi

## EUR_subnonEUR
# sample size
if [[ "${subpop}" == "EUR" ]]; then
sample_size1=233700
sample_size2=100000
else
sample_size1=311600
sample_size2=75000
fi

module load miniconda
conda activate sdprx_env

cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/

if [[ "${subpop}" != "EUR" ]]; then
# EUR_subnonEUR
# popcorn
file=/gpfs/gibbs/pi/zhao/xz674/result/sim_result/IndPRS/popcorn/sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${pop1}_sub${pop2}_${sample1_popcorn}_${sample2}_popcorn_real_corr_fold${ff}.txt
rho_est=$(grep '^pge' "${file}" | awk '{printf "%.2f", $2}')
# sdprx
EUR_sumstat="data/sim_data/summary_data/discover_validate/SDPRX/${pop1}_sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${sample1_popcorn}_SDPRX_real_all.txt"
non_EUR_sumstat="/gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/summary_data/SDPRX/${pop2}_sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${sample2}_SDPRX_real_fold${ff}.txt"
out_name="sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${pop1}_sub${pop2}_${sample1}_${sample2}_SDPRX_real_fold${ff}"
else
# subEUR_nonEUR
# popcorn
file=/gpfs/gibbs/pi/zhao/xz674/result/sim_result/IndPRS/popcorn/sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_sub${pop1}_${pop2}_${sample1_popcorn}_${sample2}_popcorn_real_corr_fold${ff}.txt
rho_est=$(grep '^pge' "${file}" | awk '{printf "%.2f", $2}')
# sdprx
EUR_sumstat="/gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/summary_data/SDPRX/${pop1}_sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${sample1_popcorn}_SDPRX_real_fold${ff}.txt"
non_EUR_sumstat="data/sim_data/summary_data/discover_validate/SDPRX/${pop2}_sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${sample2}_SDPRX_real_all.txt"
out_name="sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_sub${pop1}_${pop2}_${sample1}_${sample2}_SDPRX_real_fold${ff}"
fi

out_dir="/gpfs/gibbs/pi/zhao/xz674/result/sim_result/IndPRS/SDPRX"
mkdir -p ${out_dir}

# IFS=',' read -r -a chrs <<< "$chrs"

# for chr in "${chrs[@]}"; do

# if [ ! -f "${out_dir}/${out_name}_chr${chr}_2.txt" ] || [ $(stat -c%s "${out_dir}/${out_name}_chr${chr}_2.txt") -lt 102400 ]; then

python /gpfs/gibbs/pi/zhao/lx94/JointPRS/method/SDPRX/SDPRX.py \
--load_ld /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ref_data/SDPRX/EUR_${pop2} \
--valid /gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/data/sim_data/geno_data/All/All_test.bim \
--ss1 ${EUR_sumstat} \
--ss2 ${non_EUR_sumstat} \
--N1 ${sample_size1} \
--N2 ${sample_size2} \
--mcmc_samples 2000 \
--burn 1000 \
--force_shared True \
--chr ${chr} \
--rho ${rho_est} \
--out ${out_dir}/${out_name}_chr${chr}

# fi

# done


