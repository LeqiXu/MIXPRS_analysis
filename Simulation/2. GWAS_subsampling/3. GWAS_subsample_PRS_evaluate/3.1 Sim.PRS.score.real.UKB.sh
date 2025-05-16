# JointPRS score
job_file="/gpfs/gibbs/pi/zhao/lx94/SWIFT/code/simulation/PRS/JointPRS_score.txt"
> $job_file  # Empty the job file if it already exists

h2=0.4
rhog=0.8

sample1=ukbb
sample2=100K

for sim_i in {1..5}; do
for p in 0.001 0.01 5e-04 0.1; do
for pop2 in EUR EAS AFR SAS AMR; do

if [[ ! -e "/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/JointPRS/sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${sample1}_${sample2}_JointPRS_real_EUR_EAS_AFR_SAS_AMR_beta_${pop2}_prs_${pop2}.sscore" ]]; then

echo "module load PLINK/2; plink2 --bfile /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/${pop2}/test/${pop2} --double-id --threads 1 --score /gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/JointPRS/sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${sample1}_${sample2}_JointPRS_real_EUR_EAS_AFR_SAS_AMR_beta_${pop2}.txt header-read --out /gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/JointPRS/sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${sample1}_${sample2}_JointPRS_real_EUR_EAS_AFR_SAS_AMR_beta_${pop2}_prs_${pop2}" >> $job_file

fi

done
done
done

module load dSQ
cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/code/simulation/PRS/

dsq --job-file /gpfs/gibbs/pi/zhao/lx94/SWIFT/code/simulation/PRS/JointPRS_score.txt --partition=scavenge,day --requeue --mem=10G --cpus-per-task=1 --ntasks=1 --nodes=1 --time=24:00:00 --mail-type=ALL
sbatch dsq-JointPRS_score-$(date +%Y-%m-%d).sh


# SDPRX score
job_file="/gpfs/gibbs/pi/zhao/lx94/SWIFT/code/simulation/PRS/SDPRX_score.txt"
> $job_file  # Empty the job file if it already exists

h2=0.4
rhog=0.8

pop1=EUR
sample1=ukbb
sample2=100K

for sim_i in {1..5}; do
for p in 0.001 0.01 5e-04 0.1; do
for pop2 in EAS AFR SAS AMR; do

if [[ ! -e "/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/SDPRX/sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${sample1}_${sample2}_SDPRX_real_${pop1}_${pop2}_beta_${pop1}_prs_${pop1}.sscore" ]]; then

echo "module load PLINK/2; plink2 --bfile /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/${pop1}/test/${pop1} --double-id --threads 1 --score /gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/SDPRX/sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${sample1}_${sample2}_SDPRX_real_${pop1}_${pop2}_beta_${pop1}.txt header-read --out /gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/SDPRX/sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${sample1}_${sample2}_SDPRX_real_${pop1}_${pop2}_beta_${pop1}_prs_${pop1}" >> $job_file

fi

if [[ ! -e "/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/SDPRX/sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${sample1}_${sample2}_SDPRX_real_${pop1}_${pop2}_beta_${pop2}_prs${pop2}.sscore" ]]; then

echo "module load PLINK/2; plink2 --bfile /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/${pop2}/test/${pop2} --double-id --threads 1 --score /gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/SDPRX/sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${sample1}_${sample2}_SDPRX_real_${pop1}_${pop2}_beta_${pop2}.txt header-read --out /gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/SDPRX/sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${sample1}_${sample2}_SDPRX_real_${pop1}_${pop2}_beta_${pop2}_prs_${pop2}" >> $job_file

fi

done
done
done

module load dSQ
cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/code/simulation/PRS/

dsq --job-file /gpfs/gibbs/pi/zhao/lx94/SWIFT/code/simulation/PRS/SDPRX_score.txt --partition=scavenge,day --requeue --mem=10G --cpus-per-task=1 --ntasks=1 --nodes=1 --time=24:00:00 --mail-type=ALL
sbatch dsq-SDPRX_score-$(date +%Y-%m-%d).sh



# MIXPRS score
job_file="/gpfs/gibbs/pi/zhao/lx94/SWIFT/code/simulation/PRS/MIXPRS_score.txt"
> $job_file  # Empty the job file if it already exists

h2=0.4
rhog=0.8

pop1=EUR
sample1=ukbb
sample2=100K

for sim_i in {1..5}; do
for p in 0.001 0.01 5e-04 0.1; do
for pop2 in EUR EAS AFR SAS AMR; do

for type in full_snplist prune_snplist_1; do

if [[ "$type" == "full_snplist" ]]; then
    approx_list="FALSE"; type_original="all"
elif [[ "$type" == "prune_snplist_1" ]]; then
    approx_list="FALSE TRUE"; type_original="prune"
fi

for approx in ${approx_list}; do

if [[ ! -e "/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/MIXPRS/sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_EUR_EAS_AFR_SAS_AMR_${sample1}_${sample2}_JointPRS_SDPRX_${type}_linear_weights_approx${approx}_MIXPRS_beta_${pop2}_prs_${pop2}.sscore" ]]; then

echo "module load PLINK/2; plink2 --bfile /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/${pop2}/test/${pop2} --double-id --threads 1 --score /gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/MIXPRS/sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_EUR_EAS_AFR_SAS_AMR_${sample1}_${sample2}_JointPRS_SDPRX_${type}_linear_weights_approx${approx}_${pop2}_MIXPRS.txt header-read --out /gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/MIXPRS/sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_EUR_EAS_AFR_SAS_AMR_${sample1}_${sample2}_JointPRS_SDPRX_${type}_linear_weights_approx${approx}_MIXPRS_beta_${pop2}_prs_${pop2}" >> $job_file

fi

done

done

done
done
done

module load dSQ
cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/code/simulation/PRS/

dsq --job-file /gpfs/gibbs/pi/zhao/lx94/SWIFT/code/simulation/PRS/MIXPRS_score.txt --partition=scavenge,day --requeue --mem=10G --cpus-per-task=1 --ntasks=1 --nodes=1 --time=24:00:00 --mail-type=ALL
sbatch dsq-MIXPRS_score-$(date +%Y-%m-%d).sh
