## GLGC
GWAS_type=subsample_prune
i=1

for trait in HDL LDL TC logTG; do
for pop2 in EAS AFR SAS AMR; do
for approx in TRUE FALSE; do

if [[ ! -e "/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/MIXPRS/UKB_${trait}_MIXPRS_${GWAS_type}_${i}_linear_weights_approx${approx}_prs_${pop2}.sscore" ]]; then
        sbatch <<EOT
#!/bin/bash
#SBATCH --partition=scavenge,day
#SBATCH --requeue
#SBATCH --mem=10G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1 --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --job-name=PRS_${trait}_MIXPRS_${GWAS_type}_${i}_linear_weights_approx${approx}_prs_${pop2}
#SBATCH --output=out_PRS_${trait}_MIXPRS_${GWAS_type}_${i}_linear_weights_approx${approx}_prs_${pop2}.txt

module load PLINK/2

cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/MIXPRS/

plink2 --bfile /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/geno_data/${pop2} \
--double-id \
--threads 1 \
--score /gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/Final_weight/no_val/${trait}_${GWAS_type}_${i}_linear_weights_approx${approx}_${pop2}_MIXPRS.txt \
--out UKB_${trait}_MIXPRS_${GWAS_type}_${i}_linear_weights_approx${approx}_prs_${pop2}

EOT
fi
done
done
done

## PAGE
GWAS_type=subsample_prune
i=1

for trait in Height BMI SBP DBP PLT; do
for pop2 in EAS AFR; do
for approx in TRUE FALSE; do

if [[ ! -e "/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/MIXPRS/UKB_${trait}_MIXPRS_${GWAS_type}_${i}_linear_weights_approx${approx}_prs_${pop2}.sscore" ]]; then
        sbatch <<EOT
#!/bin/bash
#SBATCH --partition=scavenge,day
#SBATCH --requeue
#SBATCH --mem=10G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1 --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --job-name=PRS_${trait}_MIXPRS_${GWAS_type}_${i}_linear_weights_approx${approx}_prs_${pop2}
#SBATCH --output=out_PRS_${trait}_MIXPRS_${GWAS_type}_${i}_linear_weights_approx${approx}_prs_${pop2}.txt

module load PLINK/2

cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/MIXPRS/

plink2 --bfile /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/geno_data/${pop2} \
--double-id \
--threads 1 \
--score /gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/Final_weight/no_val/${trait}_${GWAS_type}_${i}_linear_weights_approx${approx}_${pop2}_MIXPRS.txt \
--out UKB_${trait}_MIXPRS_${GWAS_type}_${i}_linear_weights_approx${approx}_prs_${pop2}

EOT
fi
done
done
done

## BBJ
GWAS_type=subsample_prune
i=1

for trait in WBC NEU LYM MON EOS RBC HCT MCH MCV HB ALT ALP GGT; do
for pop2 in EAS; do
for approx in TRUE FALSE; do

if [[ ! -e "/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/MIXPRS/UKB_${trait}_MIXPRS_${GWAS_type}_${i}_linear_weights_approx${approx}_prs_${pop2}.sscore" ]]; then
        sbatch <<EOT
#!/bin/bash
#SBATCH --partition=scavenge,day
#SBATCH --requeue
#SBATCH --mem=10G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1 --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --job-name=PRS_${trait}_MIXPRS_${GWAS_type}_${i}_linear_weights_approx${approx}_prs_${pop2}
#SBATCH --output=out_PRS_${trait}_MIXPRS_${GWAS_type}_${i}_linear_weights_approx${approx}_prs_${pop2}.txt

module load PLINK/2

cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/MIXPRS/

plink2 --bfile /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/geno_data/${pop2} \
--double-id \
--threads 1 \
--score /gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/Final_weight/no_val/${trait}_${GWAS_type}_${i}_linear_weights_approx${approx}_${pop2}_MIXPRS.txt \
--out UKB_${trait}_MIXPRS_${GWAS_type}_${i}_linear_weights_approx${approx}_prs_${pop2}

EOT
fi
done
done
done

## Binary_3
GWAS_type=subsample_prune
i=1

for trait in T2D BrC; do
for pop2 in EAS AFR; do
for approx in TRUE FALSE; do

if [[ ! -e "/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/MIXPRS/UKB_${trait}_MIXPRS_${GWAS_type}_${i}_linear_weights_approx${approx}_prs_${pop2}.sscore" ]]; then
        sbatch <<EOT
#!/bin/bash
#SBATCH --partition=scavenge,day
#SBATCH --requeue
#SBATCH --mem=10G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1 --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --job-name=PRS_${trait}_MIXPRS_${GWAS_type}_${i}_linear_weights_approx${approx}_prs_${pop2}
#SBATCH --output=out_PRS_${trait}_MIXPRS_${GWAS_type}_${i}_linear_weights_approx${approx}_prs_${pop2}.txt

module load PLINK/2

cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/MIXPRS/

plink2 --bfile /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/geno_data/${pop2} \
--double-id \
--threads 1 \
--score /gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/Final_weight/no_val/${trait}_${GWAS_type}_${i}_linear_weights_approx${approx}_${pop2}_MIXPRS.txt \
--out UKB_${trait}_MIXPRS_${GWAS_type}_${i}_linear_weights_approx${approx}_prs_${pop2}

EOT
fi
done
done
done

## Binary_2
GWAS_type=subsample_prune
i=1

for trait in CAD LuC; do
for pop2 in EAS; do
for approx in TRUE FALSE; do

if [[ ! -e "/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/MIXPRS/UKB_${trait}_MIXPRS_${GWAS_type}_${i}_linear_weights_approx${approx}_prs_${pop2}.sscore" ]]; then
        sbatch <<EOT
#!/bin/bash
#SBATCH --partition=scavenge,day
#SBATCH --requeue
#SBATCH --mem=10G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1 --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --job-name=PRS_${trait}_MIXPRS_${GWAS_type}_${i}_linear_weights_approx${approx}_prs_${pop2}
#SBATCH --output=out_PRS_${trait}_MIXPRS_${GWAS_type}_${i}_linear_weights_approx${approx}_prs_${pop2}.txt

module load PLINK/2

cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/MIXPRS/

plink2 --bfile /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/geno_data/${pop2} \
--double-id \
--threads 1 \
--score /gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/Final_weight/no_val/${trait}_${GWAS_type}_${i}_linear_weights_approx${approx}_${pop2}_MIXPRS.txt \
--out UKB_${trait}_MIXPRS_${GWAS_type}_${i}_linear_weights_approx${approx}_prs_${pop2}

EOT
fi
done
done
done