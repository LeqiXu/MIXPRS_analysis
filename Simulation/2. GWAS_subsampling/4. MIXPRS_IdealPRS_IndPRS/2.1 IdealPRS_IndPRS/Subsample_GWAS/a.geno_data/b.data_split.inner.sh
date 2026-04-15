#!/bin/bash
#SBATCH --partition=scavenge,day
#SBATCH --requeue
#SBATCH --mem=30G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1 --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --job-name=data_split
#SBATCH --output=/gpfs/gibbs/pi/zhao/xz674/logs/subsample_GWAS/1022data_split/data_split_%A_%a.txt

params_file="/gpfs/gibbs/pi/zhao/xz674/simulation/a.data_prepare/subsample_GWAS/a.geno_data/b.params.txt"
IFS=' ' read pop ff <<< $(sed -n "${SLURM_ARRAY_TASK_ID}p" ${params_file})

mkdir -p /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/${pop}/discover
mkdir -p /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/${pop}/validate
mkdir -p /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/${pop}/discover_validate
mkdir -p /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/${pop}/test

# sample size
ssize="100K"

# if [[ ${pop} == "EUR" ]]; then
# ssize="100K"
# else
# ssize="20K"
# fi

## 1. Split the simulated data into discover and validate
module load PLINK/2

# ## discover set
# plink2 --bfile /gpfs/gibbs/pi/zhao/xz674/data/sim_data/geno_data/${pop}/All/${pop} \
# --double-id \
# --keep-allele-order \
# --keep /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/info/discover_${ssize}_id_fold${ff}.tsv \
# --extract /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/info/snp_rs_id_infor_hm3.snplist \
# --make-bed \
# --out /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/${pop}/discover/${pop}_${ssize}_fold${ff}

# if [[ ${pop} != "EUR" ]]; then
# plink2 --bfile /gpfs/gibbs/pi/zhao/xz674/data/sim_data/geno_data/${pop}/All/${pop} \
# --double-id \
# --keep-allele-order \
# --keep /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/info/discover_20K_id_fold${ff}.tsv \
# --extract /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/info/snp_rs_id_infor_hm3.snplist \
# --make-bed \
# --out /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/${pop}/discover/${pop}_20K_fold${ff}
# fi

# plink2 --bfile /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/1000g_phase3_data/1000G_phase3_common_norel \
# --double-id \
# --keep-allele-order \
# --keep /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/info/${pop}_discover_1kg_id_fold${ff}.tsv \
# --extract /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/1000g_phase3_data/ancestry_info/${pop}_snp.tsv \
# --make-bed \
# --out /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/${pop}/discover/${pop}_1kg_fold${ff}

# plink2 --bfile /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/1000g_phase3_data/1000G_phase3_common_norel \
# --double-id \
# --keep-allele-order \
# --keep /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/info/${pop}_discover_1kg1_id_fold${ff}.tsv \
# --extract /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/1000g_phase3_data/ancestry_info/${pop}_snp.tsv \
# --make-bed \
# --out /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/${pop}/discover/${pop}_1kg1_fold${ff}

# plink2 --bfile /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/1000g_phase3_data/1000G_phase3_common_norel \
# --double-id \
# --keep-allele-order \
# --keep /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/info/${pop}_discover_1kg2_id_fold${ff}.tsv \
# --extract /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/1000g_phase3_data/ancestry_info/${pop}_snp.tsv \
# --make-bed \
# --out /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/${pop}/discover/${pop}_1kg2_fold${ff}

## validate set
plink2 --bfile /gpfs/gibbs/pi/zhao/xz674/data/sim_data/geno_data/${pop}/All/${pop} \
--double-id \
--keep-allele-order \
--keep /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/info/validate_${ssize}_id_fold${ff}.tsv \
--extract /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/info/snp_rs_id_infor_hm3.snplist \
--make-bed \
--out /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/${pop}/validate/${pop}_${ssize}_fold${ff}

# if [[ ${pop} != "EUR" ]]; then
# plink2 --bfile /gpfs/gibbs/pi/zhao/xz674/data/sim_data/geno_data/${pop}/All/${pop} \
# --double-id \
# --keep-allele-order \
# --keep /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/info/validate_20K_id_fold${ff}.tsv \
# --extract /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/info/snp_rs_id_infor_hm3.snplist \
# --make-bed \
# --out /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/${pop}/validate/${pop}_20K_fold${ff}
# fi

# plink2 --bfile /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/1000g_phase3_data/1000G_phase3_common_norel \
# --double-id \
# --keep-allele-order \
# --keep /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/info/${pop}_validate_1kg_id_fold${ff}.tsv \
# --extract /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/1000g_phase3_data/ancestry_info/${pop}_snp.tsv \
# --make-bed \
# --out /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/${pop}/validate/${pop}_1kg_fold${ff}

# plink2 --bfile /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/1000g_phase3_data/1000G_phase3_common_norel \
# --double-id \
# --keep-allele-order \
# --keep /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/info/${pop}_validate_1kg1_id_fold${ff}.tsv \
# --extract /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/1000g_phase3_data/ancestry_info/${pop}_snp.tsv \
# --make-bed \
# --out /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/${pop}/validate/${pop}_1kg1_fold${ff}

# plink2 --bfile /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/1000g_phase3_data/1000G_phase3_common_norel \
# --double-id \
# --keep-allele-order \
# --keep /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/info/${pop}_validate_1kg2_id_fold${ff}.tsv \
# --extract /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/1000g_phase3_data/ancestry_info/${pop}_snp.tsv \
# --make-bed \
# --out /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/${pop}/validate/${pop}_1kg2_fold${ff}

# ## discover_validate set
# if [[ ${ff} == 1 ]]; then
# plink2 --bfile /gpfs/gibbs/pi/zhao/xz674/data/sim_data/geno_data/${pop}/All/${pop} \
# --double-id \
# --keep-allele-order \
# --keep /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/info/discover_validate_${ssize}_id_all.tsv \
# --extract /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/info/snp_rs_id_infor_hm3.snplist \
# --make-bed \
# --out /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/${pop}/discover_validate/${pop}_${ssize}_all

# if [[ ${pop} != "EUR" ]]; then
# plink2 --bfile /gpfs/gibbs/pi/zhao/xz674/data/sim_data/geno_data/${pop}/All/${pop} \
# --double-id \
# --keep-allele-order \
# --keep /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/info/discover_validate_20K_id_all.tsv \
# --extract /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/info/snp_rs_id_infor_hm3.snplist \
# --make-bed \
# --out /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/${pop}/discover_validate/${pop}_20K_all
# fi

# plink2 --bfile /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/1000g_phase3_data/1000G_phase3_common_norel \
# --double-id \
# --keep-allele-order \
# --keep /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/info/${pop}_discover_validate_1kg_id_all.tsv \
# --extract /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/1000g_phase3_data/ancestry_info/${pop}_snp.tsv \
# --make-bed \
# --out /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/${pop}/discover_validate/${pop}_1kg_all

# plink2 --bfile /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/1000g_phase3_data/1000G_phase3_common_norel \
# --double-id \
# --keep-allele-order \
# --keep /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/info/${pop}_discover_validate_1kg1_id_all.tsv \
# --extract /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/1000g_phase3_data/ancestry_info/${pop}_snp.tsv \
# --make-bed \
# --out /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/${pop}/discover_validate/${pop}_1kg1_all

# plink2 --bfile /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/1000g_phase3_data/1000G_phase3_common_norel \
# --double-id \
# --keep-allele-order \
# --keep /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/info/${pop}_discover_validate_1kg2_id_all.tsv \
# --extract /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/1000g_phase3_data/ancestry_info/${pop}_snp.tsv \
# --make-bed \
# --out /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/${pop}/discover_validate/${pop}_1kg2_all
# fi

# ## test set
# if [[ ${ff} == 1 ]]; then
# plink2 --bfile /gpfs/gibbs/pi/zhao/xz674/data/sim_data/geno_data/${pop}/All/${pop} \
# --double-id \
# --keep-allele-order \
# --keep /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/info/test_id.tsv \
# --extract /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/info/snp_rs_id_infor_hm3.snplist \
# --make-bed \
# --out /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/${pop}/test/${pop}
# fi


