## subsample_prune
## GLGC
GWAS_type=subsample_prune
i=1
approx=TRUE

for trait in HDL LDL TC logTG; do
for pop2 in EAS AFR SAS AMR; do
for selection_criterion in NO AIC BIC NNLS; do

if [[ "$selection_criterion" == "NO" ]]; then
    weight_name="linear_weights_approx${approx}"
elif [[ "$selection_criterion" == "AIC" || "$selection_criterion" == "BIC" ]]; then
    weight_name="linear_weights_with_forward_selection_${selection_criterion}_maxfeaturesNone_approx${approx}"
elif [[ "$selection_criterion" == "NNLS" ]]; then
    weight_name="non_negative_linear_weights_approx${approx}"
fi

if [[ ! -e "/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/MIXPRS/UKB_${trait}_MIXPRS_${GWAS_type}_${i}_${weight_name}_prs_${pop2}.sscore" ]]; then
        sbatch <<EOT
#!/bin/bash
#SBATCH --partition=scavenge,day
#SBATCH --requeue
#SBATCH --mem=10G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1 --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --job-name=PRS_${trait}_MIXPRS_${GWAS_type}_${i}_${weight_name}_prs_${pop2}
#SBATCH --output=out_PRS_${trait}_MIXPRS_${GWAS_type}_${i}_${weight_name}_prs_${pop2}.txt

module load PLINK/2

cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/MIXPRS/

plink2 --bfile /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/geno_data/${pop2} \
--double-id \
--threads 1 \
--score /gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/Final_weight/no_val/${trait}_${GWAS_type}_${i}_${weight_name}_${pop2}_MIXPRS.txt \
--out UKB_${trait}_MIXPRS_${GWAS_type}_${i}_${weight_name}_prs_${pop2}

EOT
fi
done
done
done

## PAGE
GWAS_type=subsample_prune
i=1
approx=TRUE

for trait in Height BMI SBP DBP PLT; do
for pop2 in EAS AFR; do
for selection_criterion in NO AIC BIC NNLS; do

if [[ "$selection_criterion" == "NO" ]]; then
    weight_name="linear_weights_approx${approx}"
elif [[ "$selection_criterion" == "AIC" || "$selection_criterion" == "BIC" ]]; then
    weight_name="linear_weights_with_forward_selection_${selection_criterion}_maxfeaturesNone_approx${approx}"
elif [[ "$selection_criterion" == "NNLS" ]]; then
    weight_name="non_negative_linear_weights_approx${approx}"
fi

if [[ ! -e "/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/MIXPRS/UKB_${trait}_MIXPRS_${GWAS_type}_${i}_${weight_name}_prs_${pop2}.sscore" ]]; then
        sbatch <<EOT
#!/bin/bash
#SBATCH --partition=scavenge,day
#SBATCH --requeue
#SBATCH --mem=10G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1 --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --job-name=PRS_${trait}_MIXPRS_${GWAS_type}_${i}_${weight_name}_prs_${pop2}
#SBATCH --output=out_PRS_${trait}_MIXPRS_${GWAS_type}_${i}_${weight_name}_prs_${pop2}.txt

module load PLINK/2

cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/MIXPRS/

plink2 --bfile /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/geno_data/${pop2} \
--double-id \
--threads 1 \
--score /gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/Final_weight/no_val/${trait}_${GWAS_type}_${i}_${weight_name}_${pop2}_MIXPRS.txt \
--out UKB_${trait}_MIXPRS_${GWAS_type}_${i}_${weight_name}_prs_${pop2}

EOT
fi
done
done
done

## BBJ
GWAS_type=subsample_prune
i=1
approx=TRUE

for trait in WBC NEU LYM MON EOS RBC HCT MCH MCV HB ALT ALP GGT; do
for pop2 in EAS; do
for selection_criterion in NO AIC BIC NNLS; do

if [[ "$selection_criterion" == "NO" ]]; then
    weight_name="linear_weights_approx${approx}"
elif [[ "$selection_criterion" == "AIC" || "$selection_criterion" == "BIC" ]]; then
    weight_name="linear_weights_with_forward_selection_${selection_criterion}_maxfeaturesNone_approx${approx}"
elif [[ "$selection_criterion" == "NNLS" ]]; then
    weight_name="non_negative_linear_weights_approx${approx}"
fi

if [[ ! -e "/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/MIXPRS/UKB_${trait}_MIXPRS_${GWAS_type}_${i}_${weight_name}_prs_${pop2}.sscore" ]]; then
        sbatch <<EOT
#!/bin/bash
#SBATCH --partition=scavenge,day
#SBATCH --requeue
#SBATCH --mem=10G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1 --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --job-name=PRS_${trait}_MIXPRS_${GWAS_type}_${i}_${weight_name}_prs_${pop2}
#SBATCH --output=out_PRS_${trait}_MIXPRS_${GWAS_type}_${i}_${weight_name}_prs_${pop2}.txt

module load PLINK/2

cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/MIXPRS/

plink2 --bfile /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/geno_data/${pop2} \
--double-id \
--threads 1 \
--score /gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/Final_weight/no_val/${trait}_${GWAS_type}_${i}_${weight_name}_${pop2}_MIXPRS.txt \
--out UKB_${trait}_MIXPRS_${GWAS_type}_${i}_${weight_name}_prs_${pop2}

EOT
fi
done
done
done

## Binary_3
GWAS_type=subsample_prune
i=1
approx=TRUE

for trait in T2D BrC; do
for pop2 in EAS AFR; do
for selection_criterion in NO AIC BIC NNLS; do

if [[ "$selection_criterion" == "NO" ]]; then
    weight_name="linear_weights_approx${approx}"
elif [[ "$selection_criterion" == "AIC" || "$selection_criterion" == "BIC" ]]; then
    weight_name="linear_weights_with_forward_selection_${selection_criterion}_maxfeaturesNone_approx${approx}"
elif [[ "$selection_criterion" == "NNLS" ]]; then
    weight_name="non_negative_linear_weights_approx${approx}"
fi

if [[ ! -e "/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/MIXPRS/UKB_${trait}_MIXPRS_${GWAS_type}_${i}_${weight_name}_prs_${pop2}.sscore" ]]; then
        sbatch <<EOT
#!/bin/bash
#SBATCH --partition=scavenge,day
#SBATCH --requeue
#SBATCH --mem=10G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1 --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --job-name=PRS_${trait}_MIXPRS_${GWAS_type}_${i}_${weight_name}_prs_${pop2}
#SBATCH --output=out_PRS_${trait}_MIXPRS_${GWAS_type}_${i}_${weight_name}_prs_${pop2}.txt

module load PLINK/2

cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/MIXPRS/

plink2 --bfile /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/geno_data/${pop2} \
--double-id \
--threads 1 \
--score /gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/Final_weight/no_val/${trait}_${GWAS_type}_${i}_${weight_name}_${pop2}_MIXPRS.txt \
--out UKB_${trait}_MIXPRS_${GWAS_type}_${i}_${weight_name}_prs_${pop2}

EOT
fi
done
done
done

## Binary_2
GWAS_type=subsample_prune
i=1
approx=TRUE

for trait in CAD LuC; do
for pop2 in EAS; do
for selection_criterion in NO AIC BIC NNLS; do

if [[ "$selection_criterion" == "NO" ]]; then
    weight_name="linear_weights_approx${approx}"
elif [[ "$selection_criterion" == "AIC" || "$selection_criterion" == "BIC" ]]; then
    weight_name="linear_weights_with_forward_selection_${selection_criterion}_maxfeaturesNone_approx${approx}"
elif [[ "$selection_criterion" == "NNLS" ]]; then
    weight_name="non_negative_linear_weights_approx${approx}"
fi

if [[ ! -e "/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/MIXPRS/UKB_${trait}_MIXPRS_${GWAS_type}_${i}_${weight_name}_prs_${pop2}.sscore" ]]; then
        sbatch <<EOT
#!/bin/bash
#SBATCH --partition=scavenge,day
#SBATCH --requeue
#SBATCH --mem=10G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1 --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --job-name=PRS_${trait}_MIXPRS_${GWAS_type}_${i}_${weight_name}_prs_${pop2}
#SBATCH --output=out_PRS_${trait}_MIXPRS_${GWAS_type}_${i}_${weight_name}_prs_${pop2}.txt

module load PLINK/2

cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/MIXPRS/

plink2 --bfile /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/geno_data/${pop2} \
--double-id \
--threads 1 \
--score /gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/Final_weight/no_val/${trait}_${GWAS_type}_${i}_${weight_name}_${pop2}_MIXPRS.txt \
--out UKB_${trait}_MIXPRS_${GWAS_type}_${i}_${weight_name}_prs_${pop2}

EOT
fi
done
done
done