#!/bin/bash
#SBATCH --partition=scavenge,day
#SBATCH --requeue
#SBATCH --mem=20G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1 --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --job-name=ld_calc
#SBATCH --output=/gpfs/gibbs/pi/zhao/xz674/logs/subsample_GWAS/0923ld_calc_3/ld_calc_%A_%a.txt

params_file="/gpfs/gibbs/pi/zhao/xz674/simulation/a.data_prepare/subsample_GWAS/a.geno_data/d.params.txt"
IFS=' ' read pop split geno ff <<< $(sed -n "${SLURM_ARRAY_TASK_ID}p" ${params_file})

if [[ ${ff} == "0" ]]; then
cv="all"
geno_name="${geno}"
else
cv="fold${ff}"
geno_name="${geno}_${cv}"
fi

mkdir -p /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/ref_data/ldblk/${split}/${geno_name}/${pop}

## 1. Iterate through each file in the block-wise snplist directory
module load PLINK/1.9b_6.21-x86_64

for file in "/gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/ref_data/ldblk_snplist/${split}/${geno_name}/${pop}"/*; do
# Extract the block number from the filename
blk=$(basename "${file}")

plink --bfile /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/${pop}/${split}/${pop}_${geno}_${cv} \
--keep-allele-order \
--extract ${file} \
--r square \
--out /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/ref_data/ldblk/${split}/${geno_name}/${pop}/ldblk${blk}

done
