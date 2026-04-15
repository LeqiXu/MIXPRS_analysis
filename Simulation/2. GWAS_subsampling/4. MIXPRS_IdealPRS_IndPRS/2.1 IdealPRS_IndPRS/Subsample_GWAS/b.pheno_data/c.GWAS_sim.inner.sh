#!/bin/bash
#SBATCH --partition=scavenge,day,week
#SBATCH --requeue
#SBATCH --mem=10G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1 --nodes=1
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --job-name=GWAS_simulation
#SBATCH --output=/gpfs/gibbs/pi/zhao/xz674/logs/subsample_GWAS/1004GWAS_sim/GWAS_simulation_%A_%a.txt

params_file="/gpfs/gibbs/pi/zhao/xz674/simulation/a.data_prepare/subsample_GWAS/b.pheno_data/c.params.txt"
IFS=' ' read pop i h2 p rhog ff <<< $(sed -n "${SLURM_ARRAY_TASK_ID}p" ${params_file})

mkdir -p /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/summary_data/${pop}/discover
mkdir -p /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/summary_data/${pop}/validate
mkdir -p /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/summary_data/${pop}/discover_validate

# sample size
if [[ ${pop} == "EUR" ]]; then
declare -a ssizes=("100K")
else
declare -a ssizes=("20K" "100K")
fi

## 1. Generate summary statistics
module load PLINK/1

## discover
# for ssize in "${ssizes[@]}"; do

# plink --bfile /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/${pop}/discover_validate/${pop}_100K_all \
# --extract /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/info/snp_common_hm3.snplist \
# --keep /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/info/discover_${ssize}_id_fold${ff}.tsv \
# --pheno /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/pheno_data/${pop}/discover_validate/${pop}_sim${i}_h2${h2}_p${p}_rhog${rhog}_100K_all.phen \
# --linear --allow-no-sex \
# --out /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/summary_data/${pop}/discover/${pop}_sim${i}_h2${h2}_p${p}_rhog${rhog}_${ssize}_fold${ff}

# done

# # 1kg
# plink --bfile /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/${pop}/discover_validate/${pop}_1kg_all \
# --extract /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/info/snp_common_hm3.snplist \
# --keep /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/info/${pop}_discover_1kg_id_fold${ff}.tsv \
# --pheno /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/pheno_data/${pop}/discover_validate/${pop}_sim${i}_h2${h2}_p${p}_rhog${rhog}_1kg_all.phen \
# --linear --allow-no-sex \
# --out /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/summary_data/${pop}/discover/${pop}_sim${i}_h2${h2}_p${p}_rhog${rhog}_1kg_fold${ff}

# plink --bfile /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/${pop}/discover_validate/${pop}_1kg_all \
# --extract /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/info/snp_common_hm3.snplist \
# --keep /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/info/${pop}_discover_1kg1_id_fold${ff}.tsv \
# --pheno /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/pheno_data/${pop}/discover_validate/${pop}_sim${i}_h2${h2}_p${p}_rhog${rhog}_1kg_all.phen \
# --linear --allow-no-sex \
# --out /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/summary_data/${pop}/discover/${pop}_sim${i}_h2${h2}_p${p}_rhog${rhog}_1kg1_fold${ff}

# plink --bfile /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/${pop}/discover_validate/${pop}_1kg_all \
# --extract /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/info/snp_common_hm3.snplist \
# --keep /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/info/${pop}_discover_1kg2_id_fold${ff}.tsv \
# --pheno /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/pheno_data/${pop}/discover_validate/${pop}_sim${i}_h2${h2}_p${p}_rhog${rhog}_1kg_all.phen \
# --linear --allow-no-sex \
# --out /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/summary_data/${pop}/discover/${pop}_sim${i}_h2${h2}_p${p}_rhog${rhog}_1kg2_fold${ff}

# ukbb for EUR (no cv)
# if [[ ${pop} == "EUR" ]] && [[ ${ff} == 1 ]]; then
# plink --bfile /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/geno_data/${pop} \
# --extract /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/info/snp_common_hm3.snplist \
# --keep /gpfs/gibbs/pi/zhao/xz674/data/ukbb_data/geno_data/info/discover_id.tsv \
# --pheno /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/pheno_data/${pop}/discover_validate/${pop}_sim${i}_h2${h2}_p${p}_rhog${rhog}_ukbb_all.phen \
# --linear --allow-no-sex \
# --out /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/summary_data/${pop}/discover/${pop}_sim${i}_h2${h2}_p${p}_rhog${rhog}_ukbb
# fi

# ukbb for EUR (with cv)
plink --bfile /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/geno_data/${pop} \
--extract /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/info/snp_common_hm3.snplist \
--keep /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/info/${pop}_discover_ukbb_id_fold${ff}.tsv \
--pheno /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/pheno_data/${pop}/discover_validate/${pop}_sim${i}_h2${h2}_p${p}_rhog${rhog}_ukbb_all.phen \
--linear --allow-no-sex \
--out /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/summary_data/${pop}/discover/${pop}_sim${i}_h2${h2}_p${p}_rhog${rhog}_ukbb_fold${ff}


# ## validate
# for ssize in "${ssizes[@]}"; do
# plink --bfile /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/${pop}/discover_validate/${pop}_100K_all \
# --extract /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/info/snp_common_hm3.snplist \
# --keep /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/info/validate_${ssize}_id_fold${ff}.tsv \
# --pheno /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/pheno_data/${pop}/discover_validate/${pop}_sim${i}_h2${h2}_p${p}_rhog${rhog}_100K_all.phen \
# --linear --allow-no-sex \
# --out /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/summary_data/${pop}/validate/${pop}_sim${i}_h2${h2}_p${p}_rhog${rhog}_${ssize}_fold${ff}
# done

# # 1kg
# plink --bfile /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/${pop}/discover_validate/${pop}_1kg_all \
# --extract /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/info/snp_common_hm3.snplist \
# --keep /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/info/${pop}_validate_1kg_id_fold${ff}.tsv \
# --pheno /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/pheno_data/${pop}/discover_validate/${pop}_sim${i}_h2${h2}_p${p}_rhog${rhog}_1kg_all.phen \
# --linear --allow-no-sex \
# --out /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/summary_data/${pop}/validate/${pop}_sim${i}_h2${h2}_p${p}_rhog${rhog}_1kg_fold${ff}

# plink --bfile /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/${pop}/discover_validate/${pop}_1kg_all \
# --extract /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/info/snp_common_hm3.snplist \
# --keep /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/info/${pop}_validate_1kg1_id_fold${ff}.tsv \
# --pheno /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/pheno_data/${pop}/discover_validate/${pop}_sim${i}_h2${h2}_p${p}_rhog${rhog}_1kg_all.phen \
# --linear --allow-no-sex \
# --out /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/summary_data/${pop}/validate/${pop}_sim${i}_h2${h2}_p${p}_rhog${rhog}_1kg1_fold${ff}

# plink --bfile /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/${pop}/discover_validate/${pop}_1kg_all \
# --extract /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/info/snp_common_hm3.snplist \
# --keep /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/info/${pop}_validate_1kg2_id_fold${ff}.tsv \
# --pheno /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/pheno_data/${pop}/discover_validate/${pop}_sim${i}_h2${h2}_p${p}_rhog${rhog}_1kg_all.phen \
# --linear --allow-no-sex \
# --out /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/summary_data/${pop}/validate/${pop}_sim${i}_h2${h2}_p${p}_rhog${rhog}_1kg2_fold${ff}

# # ukbb for EUR
# if [[ ${pop} == "EUR" ]] && [[ ${ff} == 1 ]]; then
# plink --bfile /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/geno_data/${pop} \
# --extract /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/info/snp_common_hm3.snplist \
# --keep /gpfs/gibbs/pi/zhao/xz674/data/ukbb_data/geno_data/info/validate_id.tsv \
# --pheno /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/pheno_data/${pop}/discover_validate/${pop}_sim${i}_h2${h2}_p${p}_rhog${rhog}_ukbb_all.phen \
# --linear --allow-no-sex \
# --out /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/summary_data/${pop}/validate/${pop}_sim${i}_h2${h2}_p${p}_rhog${rhog}_ukbb
# fi

# ## discover_validate
# if [[ ${ff} == 1 ]]; then
# for ssize in "${ssizes[@]}"; do
# plink --bfile /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/${pop}/discover_validate/${pop}_100K_all \
# --extract /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/info/snp_common_hm3.snplist \
# --keep /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/info/discover_validate_${ssize}_id_all.tsv \
# --pheno /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/pheno_data/${pop}/discover_validate/${pop}_sim${i}_h2${h2}_p${p}_rhog${rhog}_100K_all.phen \
# --linear --allow-no-sex \
# --out /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/summary_data/${pop}/discover_validate/${pop}_sim${i}_h2${h2}_p${p}_rhog${rhog}_${ssize}_all
# done

# # 1kg
# plink --bfile /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/${pop}/discover_validate/${pop}_1kg_all \
# --extract /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/info/snp_common_hm3.snplist \
# --keep /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/info/${pop}_discover_validate_1kg_id_all.tsv \
# --pheno /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/pheno_data/${pop}/discover_validate/${pop}_sim${i}_h2${h2}_p${p}_rhog${rhog}_1kg_all.phen \
# --linear --allow-no-sex \
# --out /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/summary_data/${pop}/discover_validate/${pop}_sim${i}_h2${h2}_p${p}_rhog${rhog}_1kg_all

# plink --bfile /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/${pop}/discover_validate/${pop}_1kg_all \
# --extract /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/info/snp_common_hm3.snplist \
# --keep /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/info/${pop}_discover_validate_1kg1_id_all.tsv \
# --pheno /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/pheno_data/${pop}/discover_validate/${pop}_sim${i}_h2${h2}_p${p}_rhog${rhog}_1kg_all.phen \
# --linear --allow-no-sex \
# --out /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/summary_data/${pop}/discover_validate/${pop}_sim${i}_h2${h2}_p${p}_rhog${rhog}_1kg1_all

# plink --bfile /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/${pop}/discover_validate/${pop}_1kg_all \
# --extract /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/info/snp_common_hm3.snplist \
# --keep /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/info/${pop}_discover_validate_1kg2_id_all.tsv \
# --pheno /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/pheno_data/${pop}/discover_validate/${pop}_sim${i}_h2${h2}_p${p}_rhog${rhog}_1kg_all.phen \
# --linear --allow-no-sex \
# --out /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/summary_data/${pop}/discover_validate/${pop}_sim${i}_h2${h2}_p${p}_rhog${rhog}_1kg2_all

# # ukbb for EUR
# if [[ ${pop} == "EUR" ]]; then
# plink --bfile /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/geno_data/${pop} \
# --extract /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/info/snp_common_hm3.snplist \
# --pheno /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/pheno_data/${pop}/discover_validate/${pop}_sim${i}_h2${h2}_p${p}_rhog${rhog}_ukbb_all.phen \
# --linear --allow-no-sex \
# --out /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/summary_data/${pop}/discover_validate/${pop}_sim${i}_h2${h2}_p${p}_rhog${rhog}_ukbb_all
# fi

# fi
