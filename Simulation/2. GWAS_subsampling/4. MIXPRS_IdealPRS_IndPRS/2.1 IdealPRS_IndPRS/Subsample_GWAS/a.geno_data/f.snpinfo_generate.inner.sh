#!/bin/bash
#SBATCH --partition=scavenge,day
#SBATCH --requeue
#SBATCH --mem=50G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1 --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --job-name=snpinfo_generate
#SBATCH --output=/gpfs/gibbs/pi/zhao/xz674/logs/subsample_GWAS/0923snpinfo_2/snpinfo_generate_%A_%a.txt

params_file="/gpfs/gibbs/pi/zhao/xz674/simulation/a.data_prepare/subsample_GWAS/a.geno_data/d.params.txt"
IFS=' ' read pop split geno ff <<< $(sed -n "${SLURM_ARRAY_TASK_ID}p" ${params_file})

## 1. Generate SNP info file and multi_info file
module load PLINK/1

if [[ ${geno} == "1kg" || ${geno} == "1kg1" || ${geno} == "1kg2" ]]; then
geno_type="1kg"
else
geno_type="${geno}"
fi

if [[ ${ff} == "0" ]]; then
cv="all"
else
cv="fold${ff}"
fi

# make directories
mkdir -p /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/ref_data/ldblk/${split}

# Calculate MAF using PLINK
plink --bfile /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/${pop}/${split}/${pop}_${geno}_${cv} \
--allow-extra-chr \
--freq \
--out /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/ref_data/ldblk/${split}/${pop}_${geno}_maf_${cv}

if [[ ${cv} == "all" ]]; then
# Create output file and write header
echo -e "CHR\tSNP\tBP\tA1\tA2\tMAF" > /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/ref_data/${split}/${geno}/ldblk_${geno_type}_${pop,,}/snpinfo_${geno_type}_hm3
# Merge MAF data with .bim data
awk 'NR==FNR{a1[$2]=$3; a2[$2]=$4; maf[$2]=$5; next; next} {print $1 "\t" $2 "\t" $4 "\t" a1[$2] "\t" a2[$2] "\t" maf[$2]}' /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/ref_data/ldblk/${split}/${pop}_${geno}_maf_${cv}.frq /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/${pop}/${split}/${pop}_${geno}_${cv}.bim >> /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/ref_data/${split}/${geno}/ldblk_${geno_type}_${pop,,}/snpinfo_${geno_type}_hm3
else
# Create output file and write header
echo -e "CHR\tSNP\tBP\tA1\tA2\tMAF" > /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/ref_data/${split}/${geno}_${cv}/ldblk_${geno_type}_${pop,,}/snpinfo_${geno_type}_hm3
# Merge MAF data with .bim data
awk 'NR==FNR{a1[$2]=$3; a2[$2]=$4; maf[$2]=$5; next; next} {print $1 "\t" $2 "\t" $4 "\t" a1[$2] "\t" a2[$2] "\t" maf[$2]}' /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/ref_data/ldblk/${split}/${pop}_${geno}_maf_${cv}.frq /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/${pop}/${split}/${pop}_${geno}_${cv}.bim >> /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/ref_data/${split}/${geno}_${cv}/ldblk_${geno_type}_${pop,,}/snpinfo_${geno_type}_hm3
fi

# rm /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/ref_data/ldblk/${split}/${pop}_${geno}_maf_${cv}.frq
