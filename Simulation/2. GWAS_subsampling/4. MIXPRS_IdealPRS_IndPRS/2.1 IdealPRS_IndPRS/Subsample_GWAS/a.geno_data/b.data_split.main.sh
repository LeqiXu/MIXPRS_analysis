## 0. Split dataset id
# module load R/4.2.0-foss-2020b
# Rscript /gpfs/gibbs/pi/zhao/xz674/simulation/a.data_prepare/subsample_GWAS/a.geno_data/a.id_extract.R

## 1. Data Split
# Settings
# declare -a pops=("EUR" "EAS" "AFR" "SAS" "AMR")
# declare -a ffs=({1..4})

# output_file="/gpfs/gibbs/pi/zhao/xz674/simulation/a.data_prepare/subsample_GWAS/a.geno_data/b.params.txt"
# > $output_file

# for pop in "${pops[@]}"; do
# for ff in "${ffs[@]}"; do
# echo "${pop} ${ff}" >> $output_file
# done
# done

## discover and validate job array script
params_file="/gpfs/gibbs/pi/zhao/xz674/simulation/a.data_prepare/subsample_GWAS/a.geno_data/b.params.txt"
job_count=$(wc -l < ${params_file})

sbatch --array=1-${job_count} /gpfs/gibbs/pi/zhao/xz674/simulation/a.data_prepare/subsample_GWAS/a.geno_data/b.data_split.inner.sh

## test geno data
# for pop in "${pops[@]}"; do
#     sbatch <<EOT
# #!/bin/bash
# #SBATCH --partition=scavenge,day,week
# #SBATCH --requeue
# #SBATCH --mem=30G
# #SBATCH --cpus-per-task=1
# #SBATCH --ntasks=1 --nodes=1
# #SBATCH --time=24:00:00
# #SBATCH --mail-type=ALL
# #SBATCH --job-name=data_split_${pop}
# #SBATCH --output=/gpfs/gibbs/pi/zhao/xz674/logs/subsample_GWAS/1008data_split/data_split_test_${pop}.txt

# mkdir -p /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/${pop}/test

# module load PLINK/2

# plink2 --bfile /gpfs/gibbs/pi/zhao/xz674/data/sim_data/geno_data/${pop}/All/${pop} \
# --double-id \
# --keep /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/info/test_id.tsv \
# --extract /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/info/snp_rs_id_infor_hm3.snplist \
# --make-bed \
# --out /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/${pop}/test/${pop}
# EOT
# done