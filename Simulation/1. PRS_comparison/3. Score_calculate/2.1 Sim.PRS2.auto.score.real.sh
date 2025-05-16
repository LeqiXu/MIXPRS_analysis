## 1. SDPRX
rho=0.8 

pop1=EUR
sample1=UKB

for i in {1..5}; do
for p in 0.001 0.01 5e-04 0.1; do
for sample2 in 15K 20K 25K 80K 85K 90K; do
for pop2 in EAS AFR; do
if [[ ! -e "/gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/result/sim_result/SDPRX/sim${i}_p${p}_rho${rho}/test_sim${i}_p${p}_rho${rho}_${sample1}_${sample2}_SDPRX_real_${pop1}_${pop2}_prs_${pop2}.sscore" ]]; then
if [[ "${sample2}" == "15K" || "${sample2}" == "20K" || "${sample2}" == "25K" ]]; then
sample3="15K"
elif [[ "${sample2}" == "80K" || "${sample2}" == "85K" || "${sample2}" == "90K" ]]; then
sample3="80K"
else
sample3="unknown"
fi
        sbatch <<EOT
#!/bin/bash
#SBATCH --partition=scavenge
#SBATCH --requeue
#SBATCH --mem=10G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1 --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --job-name=PRS_sim${i}_p${p}_rho${rho}_${pop1}_${pop2}_${sample1}_${sample2}_SDPRX_real
#SBATCH --output=out_PRS_sim${i}_p${p}_rho${rho}_${pop1}_${pop2}_${sample1}_${sample2}_SDPRX_real.txt

module load PLINK/2

# SDPRX
cd /gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/result/sim_result/SDPRX/sim${i}_p${p}_rho${rho}

plink2 --bfile /gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/data/sim_data/geno_data/${pop2}/test/${pop2} \
--double-id \
--threads 1 \
--extract /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/snplist_data/simulation/${pop2}_sim${i}_p${p}_rho${rho}_${sample3}_inter_snplist_real.txt \
--score /gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/result/sim_result/Final_weight/auto/SDPRX/sim${i}_p${p}_rho${rho}_${sample1}_${sample2}_SDPRX_real_${pop1}_${pop2}_beta_${pop2}.txt \
--out test_sim${i}_p${p}_rho${rho}_${sample1}_${sample2}_SDPRX_real_${pop1}_${pop2}_prs_${pop2}
EOT
fi
done
done
done
done

## 2. XPASS
rho=0.8 

pop1=EUR
sample1=UKB

for i in {1..5}; do
for p in 0.001 0.01 5e-04 0.1; do
for sample2 in 15K 20K 25K 80K 85K 90K; do
for pop2 in EAS AFR SAS AMR; do
if [[ ! -e "/gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/result/sim_result/XPASS/sim${i}_p${p}_rho${rho}/test_sim${i}_p${p}_rho${rho}_${sample1}_${sample2}_XPASS_real_${pop1}_${pop2}_prs_${pop2}.sscore" ]]; then
if [[ "${sample2}" == "15K" || "${sample2}" == "20K" || "${sample2}" == "25K" ]]; then
sample3="15K"
elif [[ "${sample2}" == "80K" || "${sample2}" == "85K" || "${sample2}" == "90K" ]]; then
sample3="80K"
else
sample3="unknown"
fi
        sbatch <<EOT
#!/bin/bash
#SBATCH --partition=scavenge
#SBATCH --requeue
#SBATCH --mem=10G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1 --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --job-name=PRS_sim${i}_p${p}_rho${rho}_${pop1}_${pop2}_${sample1}_${sample2}_XPASS_real
#SBATCH --output=out_PRS_sim${i}_p${p}_rho${rho}_${pop1}_${pop2}_${sample1}_${sample2}_XPASS_real.txt

module load PLINK/2

# XPASS
cd /gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/result/sim_result/XPASS/sim${i}_p${p}_rho${rho}

plink2 --bfile /gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/data/sim_data/geno_data/${pop2}/test/${pop2} \
--double-id \
--threads 1 \
--extract /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/snplist_data/simulation/${pop2}_sim${i}_p${p}_rho${rho}_${sample3}_inter_snplist_real.txt \
--score /gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/result/sim_result/Final_weight/auto/XPASS/sim${i}_p${p}_rho${rho}_${sample1}_${sample2}_XPASS_real_${pop1}_${pop2}_beta_${pop2}.txt \
--out test_sim${i}_p${p}_rho${rho}_${sample1}_${sample2}_XPASS_real_${pop1}_${pop2}_prs_${pop2}
EOT
fi
done
done
done
done