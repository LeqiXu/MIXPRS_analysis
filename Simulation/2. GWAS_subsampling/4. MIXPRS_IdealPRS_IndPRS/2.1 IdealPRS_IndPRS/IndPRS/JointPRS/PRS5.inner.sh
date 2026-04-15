#!/bin/bash
#SBATCH --partition=scavenge,day
#SBATCH --requeue
#SBATCH --mem=5G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1 --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --job-name=JointPRS_calc_array
#SBATCH --output=/gpfs/gibbs/pi/zhao/xz674/logs/IndPRS/1008JointPRS_calc/JointPRS_calc_%A_%a.txt

# avoid lmod error
export LMOD_IGNORE_CACHE=1

params_file="/gpfs/gibbs/pi/zhao/xz674/simulation/b.method_calculate/IndPRS/JointPRS/params.txt"
IFS=' ' read sample1 subpop sample2 h2 sim_i rhog p ff chr <<< $(sed -n "${SLURM_ARRAY_TASK_ID}p" ${params_file})

# sample size
declare -A sample_sizes
sample_sizes["EUR"]="311600"
sample_sizes["EAS"]="100000"
sample_sizes["AFR"]="100000"
sample_sizes["SAS"]="100000"
sample_sizes["AMR"]="100000"
if [[ "${subpop}" == "EUR" ]]; then
sample_sizes["${subpop}"]="233700"
else
sample_sizes["${subpop}"]="75000"
fi

out_dir="/gpfs/gibbs/pi/zhao/xz674/result/sim_result/IndPRS/JointPRS"
mkdir -p ${out_dir}

if [[ "${subpop}" == "EUR" ]]; then
out_name="sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_subEUR_EAS_AFR_SAS_AMR_${sample1}_${sample2}_JointPRS_real_fold${ff}"
elif [[ "${subpop}" == "EAS" ]]; then
out_name="sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_EUR_subEAS_AFR_SAS_AMR_${sample1}_${sample2}_JointPRS_real_fold${ff}"
elif [[ "${subpop}" == "AFR" ]]; then
out_name="sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_EUR_EAS_subAFR_SAS_AMR_${sample1}_${sample2}_JointPRS_real_fold${ff}"
elif [[ "${subpop}" == "SAS" ]]; then
out_name="sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_EUR_EAS_AFR_subSAS_AMR_${sample1}_${sample2}_JointPRS_real_fold${ff}"
elif [[ "${subpop}" == "AMR" ]]; then
out_name="sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_EUR_EAS_AFR_SAS_subAMR_${sample1}_${sample2}_JointPRS_real_fold${ff}"
fi

out_file="${out_dir}/${out_name}_AMR_pst_eff_a1_b0.5_phiauto_chr${chr}.txt"

for pop in EUR EAS AFR SAS AMR; do
  if [[ "${pop}" == "${subpop}" ]]; then
    if [[ "${pop}" == "EUR" ]]; then
        file="/gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/summary_data/PRScsx/${pop}_sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${sample1}_PRScsx_real_fold${ff}.txt"
    else
        file="/gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/summary_data/PRScsx/${pop}_sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${sample2}_PRScsx_real_fold${ff}.txt"
    fi
  else
    if [[ "${pop}" == "EUR" ]]; then
        file="data/sim_data/summary_data/discover_validate/PRScsx/${pop}_sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${sample1}_PRScsx_real_all.txt"
    else
        file="data/sim_data/summary_data/discover_validate/PRScsx/${pop}_sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${sample2}_PRScsx_real_all.txt"
    fi
  fi
  declare "${pop}_sumstat=${file}"
done

sst_file="${EUR_sumstat},${EAS_sumstat},${AFR_sumstat},${SAS_sumstat},${AMR_sumstat}"

if [ ! -f "$out_file" ] || [ $(stat -c%s "$out_file") -lt 102400 ]; then

module load miniconda
conda activate py_env

cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/

python /gpfs/gibbs/pi/zhao/lx94/JointPRS/method/JointPRS/JointPRS.py \
--ref_dir=/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ref_data/PRScsx/1kg \
--bim_prefix=/gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/data/sim_data/geno_data/All/All_test \
--sst_file=${sst_file} \
--rho_cons=1,1,1,1,1 \
--n_gwas=${sample_sizes["EUR"]},${sample_sizes["EAS"]},${sample_sizes["AFR"]},${sample_sizes["SAS"]},${sample_sizes["AMR"]} \
--chrom=${chr} \
--pop=EUR,EAS,AFR,SAS,AMR \
--out_dir=${out_dir} \
--out_name=${out_name}

fi

