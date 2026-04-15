# IdealPRS score
job_file="/gpfs/gibbs/pi/zhao/lx94/SWIFT/code/simulation/PRS/IdealPRS_score.txt"
> $job_file  # Empty the job file if it already exists

h2=0.4
rhog=0.8

pop1=EUR
sample1=ukbb
sample2=100K

approx="FALSE"

for sim_i in {1..5}; do
for p in 0.001 0.01 5e-04 0.1; do
for pop2 in EUR EAS AFR SAS AMR; do

if [[ "${pop2}" == "EUR" ]]; then
  type_list="full_snplist_ukbb"
else
  type_list="full_snplist_100K full_snplist_1kg"
fi

for type in ${type_list}; do

if [[ ! -e "/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/IdealPRS/sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_EUR_EAS_AFR_SAS_AMR_${sample1}_${sample2}_JointPRS_SDPRX_${type}_non_negative_linear_weights_approx${approx}_MIXPRS_beta_${pop2}_prs_${pop2}.sscore" ]]; then

echo "module load PLINK/2; plink2 --bfile /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/${pop2}/test/${pop2} --double-id --threads 1 --score /gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/IdealPRS/sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_EUR_EAS_AFR_SAS_AMR_${sample1}_${sample2}_JointPRS_SDPRX_${type}_non_negative_linear_weights_approx${approx}_${pop2}_MIXPRS.txt header-read --out /gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/IdealPRS/sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_EUR_EAS_AFR_SAS_AMR_${sample1}_${sample2}_JointPRS_SDPRX_${type}_non_negative_linear_weights_approx${approx}_MIXPRS_beta_${pop2}_prs_${pop2}" >> $job_file

fi

done

done
done
done

module load dSQ
cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/code/simulation/PRS/

dsq --job-file /gpfs/gibbs/pi/zhao/lx94/SWIFT/code/simulation/PRS/IdealPRS_score.txt --partition=scavenge,day,pi_zhao,week --requeue --mem=10G --cpus-per-task=1 --ntasks=1 --nodes=1 --time=24:00:00 --mail-type=ALL
sbatch dsq-IdealPRS_score-$(date +%Y-%m-%d).sh

# IndPRS score
job_file="/gpfs/gibbs/pi/zhao/lx94/SWIFT/code/simulation/PRS/IndPRS_score.txt"
> $job_file  # Empty the job file if it already exists

h2=0.4
rhog=0.8

pop1=EUR
sample1=ukbb
sample2=100K

type=full_snplist

for sim_i in {1..5}; do
for p in 0.001 0.01 5e-04 0.1; do
for pop2 in EUR EAS AFR SAS AMR; do

for weight_type in nnls lasso ridge elasticnet; do

if [[ ! -e "/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/IndPRS/sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_EUR_EAS_AFR_SAS_AMR_${sample1}_${sample2}_JointPRS_SDPRX_${type}_${weight_type}_MIXPRS_beta_${pop2}_prs_${pop2}.sscore" ]]; then

awk '$3 != 0' /gpfs/gibbs/pi/zhao/xz674/result/sim_result/IndPRS/MIXPRS/sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_EUR_EAS_AFR_SAS_AMR_${sample1}_${sample2}_JointPRS_SDPRX_${type}_${weight_type}_weights_${pop2}_MIXPRS.txt \
  > /gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/IndPRS/sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_EUR_EAS_AFR_SAS_AMR_${sample1}_${sample2}_JointPRS_SDPRX_${type}_${weight_type}_weights_${pop2}_MIXPRS.txt

echo "module load PLINK/2; plink2 --bfile /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/${pop2}/test/${pop2} --double-id --threads 1 --score /gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/IndPRS/sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_EUR_EAS_AFR_SAS_AMR_${sample1}_${sample2}_JointPRS_SDPRX_${type}_${weight_type}_weights_${pop2}_MIXPRS.txt header-read --out /gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/IndPRS/sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_EUR_EAS_AFR_SAS_AMR_${sample1}_${sample2}_JointPRS_SDPRX_${type}_${weight_type}_MIXPRS_beta_${pop2}_prs_${pop2}" >> $job_file

fi

done

done
done
done

module load dSQ
cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/code/simulation/PRS/

dsq --job-file /gpfs/gibbs/pi/zhao/lx94/SWIFT/code/simulation/PRS/IndPRS_score.txt --partition=scavenge,day --requeue --mem=10G --cpus-per-task=1 --ntasks=1 --nodes=1 --time=24:00:00 --mail-type=ALL
sbatch dsq-IndPRS_score-$(date +%Y-%m-%d).sh
