## 1. JointPRS
GWAS_type=real

for trait in T2D BrC; do
for pop2 in EAS AFR; do

if [[ ! -e "/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/JointPRS/UKB_${trait}_JointPRS_${GWAS_type}_EUR_EAS_AFR_beta_AFR_prs_${pop2}.sscore" ]]; then
        sbatch <<EOT
#!/bin/bash
#SBATCH --partition=scavenge,day
#SBATCH --requeue
#SBATCH --mem=10G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1 --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --job-name=PRS_${trait}_JointPRS_${GWAS_type}_EUR_EAS_AFR_prs_${pop2}
#SBATCH --output=out_PRS_${trait}_JointPRS_${GWAS_type}_EUR_EAS_AFR_prs_${pop2}.txt

module load PLINK/2

cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/JointPRS/

plink2 --bfile /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/geno_data/${pop2} \
--double-id \
--threads 1 \
--score /gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/result/summary_result/Final_weight/no_val/JointPRS/${trait}_JointPRS_EUR_EAS_AFR_beta_EUR.txt header-read \
--out UKB_${trait}_JointPRS_${GWAS_type}_EUR_EAS_AFR_beta_EUR_prs_${pop2}

plink2 --bfile /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/geno_data/${pop2} \
--double-id \
--threads 1 \
--score /gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/result/summary_result/Final_weight/no_val/JointPRS/${trait}_JointPRS_EUR_EAS_AFR_beta_EAS.txt header-read \
--out UKB_${trait}_JointPRS_${GWAS_type}_EUR_EAS_AFR_beta_EAS_prs_${pop2}

plink2 --bfile /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/geno_data/${pop2} \
--double-id \
--threads 1 \
--score /gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/result/summary_result/Final_weight/no_val/JointPRS/${trait}_JointPRS_EUR_EAS_AFR_beta_AFR.txt header-read \
--out UKB_${trait}_JointPRS_${GWAS_type}_EUR_EAS_AFR_beta_AFR_prs_${pop2}

EOT
fi
done
done

GWAS_type=real

for trait in CAD LuC; do
for pop2 in EAS; do

if [[ ! -e "/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/JointPRS/UKB_${trait}_JointPRS_${GWAS_type}_EUR_${pop2}_beta_EUR_prs_${pop2}.sscore" ]]; then
        sbatch <<EOT
#!/bin/bash
#SBATCH --partition=scavenge,day
#SBATCH --requeue
#SBATCH --mem=10G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1 --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --job-name=PRS_${trait}_JointPRS_${GWAS_type}_EUR_${pop2}_prs_${pop2}
#SBATCH --output=out_PRS_${trait}_JointPRS_${GWAS_type}_EUR_${pop2}_prs_${pop2}.txt

module load PLINK/2

cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/JointPRS/

plink2 --bfile /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/geno_data/${pop2} \
--double-id \
--threads 1 \
--score /gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/result/summary_result/Final_weight/no_val/JointPRS/${trait}_JointPRS_EUR_${pop2}_beta_${pop2}.txt header-read \
--out UKB_${trait}_JointPRS_${GWAS_type}_EUR_${pop2}_beta_${pop2}_prs_${pop2}

plink2 --bfile /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/geno_data/${pop2} \
--double-id \
--threads 1 \
--score /gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/result/summary_result/Final_weight/no_val/JointPRS/${trait}_JointPRS_EUR_${pop2}_beta_EUR.txt header-read \
--out UKB_${trait}_JointPRS_${GWAS_type}_EUR_${pop2}_beta_EUR_prs_${pop2}
EOT
fi
done
done


## 2. SDPRX
GWAS_type=real

for trait in T2D BrC; do
for pop2 in EAS AFR; do

if [[ ! -e "/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/SDPRX/UKB_${trait}_SDPRX_${GWAS_type}_EUR_${pop2}_beta_EUR_prs_${pop2}.sscore" ]]; then
        sbatch <<EOT
#!/bin/bash
#SBATCH --partition=scavenge,day
#SBATCH --requeue
#SBATCH --mem=10G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1 --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --job-name=PRS_${trait}_SDPRX_${GWAS_type}_EUR_${pop2}_prs_${pop2}
#SBATCH --output=out_PRS_${trait}_SDPRX_${GWAS_type}_EUR_${pop2}_prs_${pop2}.txt

module load PLINK/2

cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/SDPRX/

plink2 --bfile /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/geno_data/${pop2} \
--double-id \
--threads 1 \
--score /gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/result/summary_result/Final_weight/no_val/SDPRX/${trait}_SDPRX_EUR_${pop2}_beta_EUR.txt header-read \
--out UKB_${trait}_SDPRX_${GWAS_type}_EUR_${pop2}_beta_EUR_prs_${pop2}

plink2 --bfile /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/geno_data/${pop2} \
--double-id \
--threads 1 \
--score /gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/result/summary_result/Final_weight/no_val/SDPRX/${trait}_SDPRX_EUR_EAS_beta_EAS.txt header-read \
--out UKB_${trait}_SDPRX_${GWAS_type}_EUR_EAS_beta_EAS_prs_${pop2}

plink2 --bfile /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/geno_data/${pop2} \
--double-id \
--threads 1 \
--score /gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/result/summary_result/Final_weight/no_val/SDPRX/${trait}_SDPRX_EUR_AFR_beta_AFR.txt header-read \
--out UKB_${trait}_SDPRX_${GWAS_type}_EUR_AFR_beta_AFR_prs_${pop2}

EOT
fi
done
done

GWAS_type=real

for trait in CAD LuC; do
for pop2 in EAS; do

if [[ ! -e "/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/SDPRX/UKB_${trait}_SDPRX_${GWAS_type}_EUR_${pop2}_beta_EUR_prs_${pop2}.sscore" ]]; then
        sbatch <<EOT
#!/bin/bash
#SBATCH --partition=scavenge,day
#SBATCH --requeue
#SBATCH --mem=10G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1 --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --job-name=PRS_${trait}_SDPRX_${GWAS_type}_EUR_${pop2}_prs_${pop2}
#SBATCH --output=out_PRS_${trait}_SDPRX_${GWAS_type}_EUR_${pop2}_prs_${pop2}.txt

module load PLINK/2

cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/SDPRX/

plink2 --bfile /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/geno_data/${pop2} \
--double-id \
--threads 1 \
--score /gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/result/summary_result/Final_weight/no_val/SDPRX/${trait}_SDPRX_EUR_${pop2}_beta_${pop2}.txt header-read \
--out UKB_${trait}_SDPRX_${GWAS_type}_EUR_${pop2}_beta_${pop2}_prs_${pop2}

plink2 --bfile /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/geno_data/${pop2} \
--double-id \
--threads 1 \
--score /gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/result/summary_result/Final_weight/no_val/SDPRX/${trait}_SDPRX_EUR_${pop2}_beta_EUR.txt header-read \
--out UKB_${trait}_SDPRX_${GWAS_type}_EUR_${pop2}_beta_EUR_prs_${pop2}
EOT
fi
done
done


## 3. MIXPRS
GWAS_type=subsample_prune

for trait in T2D BrC; do
for pop2 in EAS AFR; do

if [[ ! -e "/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/MIXPRS/UKB_${trait}_MIXPRS_${GWAS_type}_prs_${pop2}.sscore" ]]; then
        sbatch <<EOT
#!/bin/bash
#SBATCH --partition=scavenge,day
#SBATCH --requeue
#SBATCH --mem=10G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1 --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --job-name=PRS_${trait}_MIXPRS_${GWAS_type}_prs_${pop2}
#SBATCH --output=out_PRS_${trait}_MIXPRS_${GWAS_type}_prs_${pop2}.txt

module load PLINK/2

cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/MIXPRS/

plink2 --bfile /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/geno_data/${pop2} \
--double-id \
--threads 1 \
--score /gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/Final_weight/no_val/${trait}_${pop2}_MIXPRS.txt \
--out UKB_${trait}_MIXPRS_${GWAS_type}_prs_${pop2}

EOT
fi
done
done


GWAS_type=subsample_prune

for trait in CAD LuC; do
for pop2 in EAS; do

if [[ ! -e "/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/MIXPRS/UKB_${trait}_MIXPRS_${GWAS_type}_prs_${pop2}.sscore" ]]; then
        sbatch <<EOT
#!/bin/bash
#SBATCH --partition=scavenge,day
#SBATCH --requeue
#SBATCH --mem=10G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1 --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --job-name=PRS_${trait}_MIXPRS_${GWAS_type}_prs_${pop2}
#SBATCH --output=out_PRS_${trait}_MIXPRS_${GWAS_type}_prs_${pop2}.txt

module load PLINK/2

cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/MIXPRS/

plink2 --bfile /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/geno_data/${pop2} \
--double-id \
--threads 1 \
--score /gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/Final_weight/no_val/${trait}_${pop2}_MIXPRS.txt \
--out UKB_${trait}_MIXPRS_${GWAS_type}_prs_${pop2}

EOT
fi
done
done
