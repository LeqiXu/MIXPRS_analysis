## 3. MIXPRS update
h2=0.4
rhog=0.8 

pop1=EUR
sample1=UKB

type=prune_snplist_1
approx=TRUE

for sim_i in {1..5}; do
for p in 0.001 0.01 5e-04 0.1; do
for sample2 in 15K 80K; do
for pop2 in EAS AFR SAS AMR; do
if [[ ! -e "/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/MIXPRS/sim${sim_i}_p${p}_rho${rhog}_${sample2}_MIXPRS_${type}_prs_${pop2}.sscore" ]]; then

        sbatch <<EOT
#!/bin/bash
#SBATCH --partition=scavenge,day
#SBATCH --requeue
#SBATCH --mem=10G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1 --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --job-name=PRS_sim${sim_i}_p${p}_rho${rhog}_EUR_EAS_AFR_SAS_AMR_${sample1}_${sample2}_MIXPRS_real
#SBATCH --output=out_PRS_sim${sim_i}_p${p}_rho${rhog}_EUR_EAS_AFR_SAS_AMR_${sample1}_${sample2}_MIXPRS_real.txt

module load PLINK/2

# MIXPRS
cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/MIXPRS

plink2 --bfile /gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/data/sim_data/geno_data/${pop2}/test/${pop2} \
--double-id \
--threads 1 \
--extract /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/snplist_data/simulation/${pop2}_sim${sim_i}_p${p}_rho${rhog}_${sample2}_inter_snplist_real.txt \
--score /gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/MIXPRS/sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_EUR_EAS_AFR_SAS_AMR_${sample1}_${sample2}_JointPRS_SDPRX_${type}_non_negative_linear_weights_approx${approx}_${pop2}_MIXPRS.txt \
--out sim${sim_i}_p${p}_rho${rhog}_${sample2}_MIXPRS_${type}_prs_${pop2}
EOT
fi
done
done
done
done