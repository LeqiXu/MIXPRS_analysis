# MIXPRS score
job_file="/gpfs/gibbs/pi/zhao/lx94/SWIFT/code/simulation/PRS/MIXPRS_score.txt"
> $job_file  # Empty the job file if it already exists

h2=0.4
rhog=0.8

pop1=EUR
sample1=ukbb
sample2=100K

type=prune_snplist_1
type_original="prune"
approx_list="TRUE"

for sim_i in {1..5}; do
for p in 0.001 0.01 5e-04 0.1; do
for pop2 in EUR EAS AFR SAS AMR; do

for approx in ${approx_list}; do

if [[ ! -e "/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/MIXPRS/sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_EUR_EAS_AFR_SAS_AMR_${sample1}_${sample2}_JointPRS_SDPRX_${type}_non_negative_linear_weights_approx${approx}_MIXPRS_beta_${pop2}_prs_${pop2}.sscore" ]]; then

echo "module load PLINK/2; plink2 --bfile /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/${pop2}/test/${pop2} --double-id --threads 1 --score /gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/MIXPRS/sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_EUR_EAS_AFR_SAS_AMR_${sample1}_${sample2}_JointPRS_SDPRX_${type}_non_negative_linear_weights_approx${approx}_${pop2}_MIXPRS.txt header-read --out /gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/MIXPRS/sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_EUR_EAS_AFR_SAS_AMR_${sample1}_${sample2}_JointPRS_SDPRX_${type}_non_negative_linear_weights_approx${approx}_MIXPRS_beta_${pop2}_prs_${pop2}" >> $job_file

fi

done

done
done
done

module load dSQ
cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/code/simulation/PRS/

dsq --job-file /gpfs/gibbs/pi/zhao/lx94/SWIFT/code/simulation/PRS/MIXPRS_score.txt --partition=scavenge,day --requeue --mem=10G --cpus-per-task=1 --ntasks=1 --nodes=1 --time=24:00:00 --mail-type=ALL
sbatch dsq-MIXPRS_score-$(date +%Y-%m-%d).sh
