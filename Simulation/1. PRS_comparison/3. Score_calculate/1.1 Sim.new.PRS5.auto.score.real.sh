## 1. JointPRS
h2=0.4
rhog=0.8 

pop1=EUR
sample1=UKB

for sim_i in {1..5}; do
for p in 0.001 0.01 5e-04 0.1; do
for sample2 in 25K 90K; do
for pop2 in EAS AFR SAS AMR; do
if [[ ! -e "/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/JointPRS/test_sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${sample1}_${sample2}_JointPRS_real_EUR_EAS_AFR_SAS_AMR_beta_${pop2}_prs_${pop2}.sscore" ]]; then
if [[ "${sample2}" == "15K" || "${sample2}" == "20K" || "${sample2}" == "25K" ]]; then
sample3="15K"
elif [[ "${sample2}" == "80K" || "${sample2}" == "85K" || "${sample2}" == "90K" ]]; then
sample3="80K"
else
sample3="unknown"
fi
        sbatch <<EOT
#!/bin/bash
#SBATCH --partition=scavenge,day,week
#SBATCH --requeue
#SBATCH --mem=10G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1 --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --job-name=PRS_sim${sim_i}_p${p}_rho${rhog}_EUR_EAS_AFR_SAS_AMR_${sample1}_${sample2}_JointPRS_real
#SBATCH --output=out_PRS_sim${sim_i}_p${p}_rho${rhog}_EUR_EAS_AFR_SAS_AMR_${sample1}_${sample2}_JointPRS_real.txt

module load PLINK/2

# EUR
cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/JointPRS

plink2 --bfile /gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/data/sim_data/geno_data/${pop2}/test/${pop2} \
--double-id \
--threads 1 \
--extract /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/snplist_data/simulation/${pop2}_sim${sim_i}_p${p}_rho${rhog}_${sample3}_inter_snplist_real.txt \
--score /gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/JointPRS/sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${sample1}_${sample2}_JointPRS_real_EUR_EAS_AFR_SAS_AMR_beta_EUR.txt \
--out test_sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${sample1}_${sample2}_JointPRS_real_EUR_EAS_AFR_SAS_AMR_beta_EUR_prs_${pop2}

# EAS
cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/JointPRS

plink2 --bfile /gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/data/sim_data/geno_data/${pop2}/test/${pop2} \
--double-id \
--threads 1 \
--extract /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/snplist_data/simulation/${pop2}_sim${sim_i}_p${p}_rho${rhog}_${sample3}_inter_snplist_real.txt \
--score /gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/JointPRS/sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${sample1}_${sample2}_JointPRS_real_EUR_EAS_AFR_SAS_AMR_beta_EAS.txt \
--out test_sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${sample1}_${sample2}_JointPRS_real_EUR_EAS_AFR_SAS_AMR_beta_EAS_prs_${pop2}

# AFR
cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/JointPRS

plink2 --bfile /gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/data/sim_data/geno_data/${pop2}/test/${pop2} \
--double-id \
--threads 1 \
--extract /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/snplist_data/simulation/${pop2}_sim${sim_i}_p${p}_rho${rhog}_${sample3}_inter_snplist_real.txt \
--score /gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/JointPRS/sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${sample1}_${sample2}_JointPRS_real_EUR_EAS_AFR_SAS_AMR_beta_AFR.txt \
--out test_sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${sample1}_${sample2}_JointPRS_real_EUR_EAS_AFR_SAS_AMR_beta_AFR_prs_${pop2}

# SAS
cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/JointPRS

plink2 --bfile /gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/data/sim_data/geno_data/${pop2}/test/${pop2} \
--double-id \
--threads 1 \
--extract /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/snplist_data/simulation/${pop2}_sim${sim_i}_p${p}_rho${rhog}_${sample3}_inter_snplist_real.txt \
--score /gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/JointPRS/sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${sample1}_${sample2}_JointPRS_real_EUR_EAS_AFR_SAS_AMR_beta_SAS.txt \
--out test_sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${sample1}_${sample2}_JointPRS_real_EUR_EAS_AFR_SAS_AMR_beta_SAS_prs_${pop2}

# AMR
cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/JointPRS

plink2 --bfile /gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/data/sim_data/geno_data/${pop2}/test/${pop2} \
--double-id \
--threads 1 \
--extract /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/snplist_data/simulation/${pop2}_sim${sim_i}_p${p}_rho${rhog}_${sample3}_inter_snplist_real.txt \
--score /gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/JointPRS/sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${sample1}_${sample2}_JointPRS_real_EUR_EAS_AFR_SAS_AMR_beta_AMR.txt \
--out test_sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${sample1}_${sample2}_JointPRS_real_EUR_EAS_AFR_SAS_AMR_beta_AMR_prs_${pop2}
EOT
fi
done
done
done
done


## 2. MIXPRS
h2=0.4
rhog=0.8 

pop1=EUR
sample1=UKB

type=prune_snplist_1
approx=TRUE

for sim_i in {1..5}; do
for p in 0.001 0.01 5e-04 0.1; do
for sample2 in 25K 90K; do
for pop2 in EAS AFR SAS AMR; do
if [[ ! -e "/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/MIXPRS/test_sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_EUR_EAS_AFR_SAS_AMR_${sample1}_${sample2}_JointPRS_SDPRX_${type}_non_negative_linear_weights_approx${approx}_${pop2}_MIXPRS_prs_${pop2}.sscore" ]]; then
if [[ "${sample2}" == "15K" || "${sample2}" == "20K" || "${sample2}" == "25K" ]]; then
sample3="15K"
elif [[ "${sample2}" == "80K" || "${sample2}" == "85K" || "${sample2}" == "90K" ]]; then
sample3="80K"
else
sample3="unknown"
fi
        sbatch <<EOT
#!/bin/bash
#SBATCH --partition=scavenge,day,week
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
--extract /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/snplist_data/simulation/${pop2}_sim${sim_i}_p${p}_rho${rhog}_${sample3}_inter_snplist_real.txt \
--score /gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/MIXPRS/sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_EUR_EAS_AFR_SAS_AMR_${sample1}_${sample2}_JointPRS_SDPRX_${type}_non_negative_linear_weights_approx${approx}_${pop2}_MIXPRS.txt \
--out test_sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_EUR_EAS_AFR_SAS_AMR_${sample1}_${sample2}_JointPRS_SDPRX_${type}_non_negative_linear_weights_approx${approx}_${pop2}_MIXPRS_prs_${pop2}
EOT
fi
done
done
done
done