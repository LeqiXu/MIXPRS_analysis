#!/bin/bash
#SBATCH --partition=scavenge,day
#SBATCH --requeue
#SBATCH --mem=3G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1 --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --job-name=popcorn_calc_array
#SBATCH --output=/gpfs/gibbs/pi/zhao/xz674/logs/IndPRS/1006popcorn_calc/popcorn_calc_%A_%a.txt

# avoid lmod error
export LMOD_IGNORE_CACHE=1

params_file="/gpfs/gibbs/pi/zhao/xz674/simulation/b.method_calculate/IndPRS/popcorn/params2.txt"
IFS=' ' read pop1 sample1 pop2 sample2 subpop h2 sim_i rhog p <<< $(sed -n "${SLURM_ARRAY_TASK_ID}p" ${params_file})

module load miniconda
conda activate popcorn

cd /gpfs/gibbs/pi/zhao/xz674/Popcorn

out_dir="/gpfs/gibbs/pi/zhao/xz674/result/sim_result/IndPRS/popcorn"
mkdir -p ${out_dir}

declare -a ffs=({1..4})

for ff in "${ffs[@]}"; do

if [[ ${subpop} == "EUR" ]]; then
sfile1="/gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/summary_data/popcorn/${pop1}_sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${sample1}_popcorn_real_fold${ff}.txt"
sfile2="/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/sim_data/summary_data/discover_validate/popcorn/${pop2}_sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${sample2}_popcorn_real_all.txt"
popcorn fit -v 0 \
--cfile /gpfs/gibbs/pi/zhao/lx94/JointPRS/method/Popcorn/ref/${pop1}_${pop2}_all_gen_eff.cscore \
--gen_effect \
--sfile1 ${sfile1} \
--sfile2 ${sfile2} \
${out_dir}/sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_sub${pop1}_${pop2}_${sample1}_${sample2}_popcorn_real_corr_fold${ff}.txt

else
sfile1="/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/sim_data/summary_data/discover_validate/popcorn/${pop1}_sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${sample1}_popcorn_real_all.txt"
sfile2="/gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/summary_data/popcorn/${pop2}_sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${sample2}_popcorn_real_fold${ff}.txt"
popcorn fit -v 0 \
--cfile /gpfs/gibbs/pi/zhao/lx94/JointPRS/method/Popcorn/ref/${pop1}_${pop2}_all_gen_eff.cscore \
--gen_effect \
--sfile1 ${sfile1} \
--sfile2 ${sfile2} \
${out_dir}/sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${pop1}_sub${pop2}_${sample1}_${sample2}_popcorn_real_corr_fold${ff}.txt

fi

done



