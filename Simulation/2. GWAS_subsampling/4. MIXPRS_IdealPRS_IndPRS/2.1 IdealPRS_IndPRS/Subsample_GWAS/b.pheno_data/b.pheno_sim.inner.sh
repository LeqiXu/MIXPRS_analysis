#!/bin/bash
#SBATCH --partition=scavenge,day
#SBATCH --requeue
#SBATCH --mem=200G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1 --nodes=1
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --job-name=phenotype_simulation
#SBATCH --output=/gpfs/gibbs/pi/zhao/xz674/logs/subsample_GWAS/0922pheno_sim_3/phenotype_simulation_%A_%a.txt

params_file="/gpfs/gibbs/pi/zhao/xz674/simulation/a.data_prepare/subsample_GWAS/b.pheno_data/b.params.txt"
IFS=' ' read pop i h2 p rhog <<< $(sed -n "${SLURM_ARRAY_TASK_ID}p" ${params_file})

mkdir -p /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/pheno_data/${pop}/discover_validate
mkdir -p /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/pheno_data/${pop}/test

# sample size
if [[ ${pop} == "EUR" ]]; then
declare -a ssizes=("100K")
else
declare -a ssizes=("20K" "100K")
fi

# discover_validate phenotype
# simulated data
for ssize in "${ssizes[@]}"; do
if [[ ! -e "/gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/pheno_data/${pop}/discover_validate/${pop}_sim${i}_h2${h2}_p${p}_rhog${rhog}_${ssize}_all.phen" ]]; then
  /gpfs/gibbs/pi/zhao/lx94/JointPRS/method/gcta64/gcta-1.94.1 --bfile /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/${pop}/discover_validate/${pop}_${ssize}_all \
  --extract /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/effect_data/${pop}_sim${i}_h2${h2}_p${p}_rhog${rhog}_simu.snplist \
  --keep /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/info/discover_validate_${ssize}_id_all.tsv \
  --simu-qt \
  --simu-causal-loci /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/effect_data/${pop}_sim${i}_h2${h2}_p${p}_rhog${rhog}_simu.txt \
  --simu-hsq ${h2} \
  --simu-rep 1 \
  --out /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/pheno_data/${pop}/discover_validate/${pop}_sim${i}_h2${h2}_p${p}_rhog${rhog}_${ssize}_all
fi
done

# 1kg data
if [[ ! -e "/gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/pheno_data/${pop}/discover_validate/${pop}_sim${i}_h2${h2}_p${p}_rhog${rhog}_1kg_all.phen" ]]; then
  /gpfs/gibbs/pi/zhao/lx94/JointPRS/method/gcta64/gcta-1.94.1 --bfile /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/${pop}/discover_validate/${pop}_1kg_all \
  --extract /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/effect_data/${pop}_sim${i}_h2${h2}_p${p}_rhog${rhog}_1kg.snplist \
  --keep /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/info/${pop}_discover_validate_1kg_id_all.tsv \
  --simu-qt \
  --simu-causal-loci /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/effect_data/${pop}_sim${i}_h2${h2}_p${p}_rhog${rhog}_1kg.txt \
  --simu-hsq ${h2} \
  --simu-rep 1 \
  --out /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/pheno_data/${pop}/discover_validate/${pop}_sim${i}_h2${h2}_p${p}_rhog${rhog}_1kg_all
fi

# ukbb data
if [[ ${pop} == "EUR" ]]; then
if [[ ! -e "/gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/pheno_data/${pop}/discover_validate/${pop}_sim${i}_h2${h2}_p${p}_rhog${rhog}_ukbb_all.phen" ]]; then
  /gpfs/gibbs/pi/zhao/lx94/JointPRS/method/gcta64/gcta-1.94.1 --bfile /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/geno_data/EUR \
  --extract /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/effect_data/${pop}_sim${i}_h2${h2}_p${p}_rhog${rhog}_ukbb.snplist \
  --simu-qt \
  --simu-causal-loci /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/effect_data/${pop}_sim${i}_h2${h2}_p${p}_rhog${rhog}_ukbb.txt \
  --simu-hsq ${h2} \
  --simu-rep 1 \
  --out /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/pheno_data/${pop}/discover_validate/${pop}_sim${i}_h2${h2}_p${p}_rhog${rhog}_ukbb_all
fi
fi

# test phenotype
if [[ ! -e "/gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/pheno_data/${pop}/test/${pop}_sim${i}_h2${h2}_p${p}_rhog${rhog}_10K.phen" ]]; then
  /gpfs/gibbs/pi/zhao/lx94/JointPRS/method/gcta64/gcta-1.94.1 --bfile /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/${pop}/test/${pop} \
  --extract /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/effect_data/${pop}_sim${i}_h2${h2}_p${p}_rhog${rhog}_simu.snplist \
  --keep /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/info/test_id.tsv \
  --simu-qt \
  --simu-causal-loci /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/effect_data/${pop}_sim${i}_h2${h2}_p${p}_rhog${rhog}_simu.txt \
  --simu-hsq ${h2} \
  --simu-rep 1 \
  --out /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/pheno_data/${pop}/test/${pop}_sim${i}_h2${h2}_p${p}_rhog${rhog}_10K
fi
