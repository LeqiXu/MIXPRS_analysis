#!/bin/bash
#SBATCH --partition=scavenge,day
#SBATCH --requeue
#SBATCH --mem=15G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1 --nodes=1
#SBATCH --time=2:00:00
#SBATCH --mail-type=ALL
#SBATCH --job-name=PRS_SDPRX_real_array
#SBATCH --output=/gpfs/gibbs/pi/zhao/xz674/logs/IndPRS/1025SDPRX_score/PRS_SDPRX_real_%A_%a.txt

# avoid lmod error
export LMOD_IGNORE_CACHE=1

params_file="/gpfs/gibbs/pi/zhao/xz674/simulation/b.method_calculate/IndPRS/linear_combine/SDPRX.params.txt"
IFS=' ' read pop1 sample1 pop2 sample2 subpop h2 sim_i rhog p target_pop <<< $(sed -n "${SLURM_ARRAY_TASK_ID}p" ${params_file})

module load PLINK/2

out_dir="/gpfs/gibbs/pi/zhao/xz674/result/sim_result/IndPRS/SDPRX"

declare -a ffs=({1..4})

for ff in "${ffs[@]}"; do

if [[ "${subpop}" != "EUR" ]]; then
# EUR_subnonEUR
out_name="sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${pop1}_sub${pop2}_${sample1}_${sample2}_SDPRX_real_fold${ff}"
else
# subEUR_nonEUR
out_name="sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_sub${pop1}_${pop2}_${sample1}_${sample2}_SDPRX_real_fold${ff}"
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
