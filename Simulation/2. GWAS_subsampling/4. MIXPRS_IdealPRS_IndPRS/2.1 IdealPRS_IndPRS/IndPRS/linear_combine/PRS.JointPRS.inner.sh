#!/bin/bash
#SBATCH --partition=scavenge,day
#SBATCH --requeue
#SBATCH --mem=15G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1 --nodes=1
#SBATCH --time=2:00:00
#SBATCH --mail-type=ALL
#SBATCH --job-name=PRS_JointPRS_real_array
#SBATCH --output=/gpfs/gibbs/pi/zhao/xz674/logs/IndPRS/1025JointPRS_score/PRS_JointPRS_real_%A_%a.txt

# avoid lmod error
export LMOD_IGNORE_CACHE=1

params_file="/gpfs/gibbs/pi/zhao/xz674/simulation/b.method_calculate/IndPRS/linear_combine/JointPRS.params.txt"
IFS=' ' read sample1 subpop sample2 h2 sim_i rhog p target_pop <<< $(sed -n "${SLURM_ARRAY_TASK_ID}p" ${params_file})

module load PLINK/2

out_dir="/gpfs/gibbs/pi/zhao/xz674/result/sim_result/IndPRS/JointPRS"

declare -a ffs=({1..4})

for ff in "${ffs[@]}"; do

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

out_file="${out_dir}/${out_name}_beta${target_pop}.txt"

if [[ ${subpop} == "EUR" ]]; then

plink2 --bfile /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/geno_data/${subpop}_hm3 \
--keep /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/info/${subpop}_validate_ukbb_id_fold${ff}.tsv \
--double-id \
--threads 1 \
--score ${out_file} \
--out ${out_dir}/${out_name}_beta_${target_pop}_prs_${target_pop}

else

plink2 --bfile /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/${subpop}/validate/${subpop}_${sample2}_fold${ff} \
--double-id \
--threads 1 \
--score ${out_file} \
--out ${out_dir}/${out_name}_beta_${target_pop}_prs_${target_pop}

fi
done
