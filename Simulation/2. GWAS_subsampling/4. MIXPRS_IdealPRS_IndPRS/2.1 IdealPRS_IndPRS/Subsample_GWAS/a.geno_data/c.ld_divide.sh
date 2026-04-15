#!/bin/bash
#SBATCH --partition=scavenge,day
#SBATCH --requeue
#SBATCH --mem=30G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1 --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --job-name=ld_divide
#SBATCH --output=/gpfs/gibbs/pi/zhao/xz674/logs/subsample_GWAS/0923ld_div_2/ld_divide.txt

## make file directory
# declare -a pops=("EUR" "EAS" "AFR" "SAS" "AMR")
# declare -a splits=("discover_validate" "discover")

# for pop in "${pops[@]}"; do
# # declare -a genos=("1kg")
# if [[ ${pop} == "EUR" ]]; then
# declare -a genos=("100K" "1kg")
# else
# declare -a genos=("100K" "20K" "1kg")
# fi
# for split in "${splits[@]}"; do
# mkdir -p /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/ref_data/ldblk_info/${split}
# for geno in "${genos[@]}"; do
# if [[ ${split} == "discover_validate" ]]; then
# mkdir -p /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/ref_data/ldblk_snplist/${split}/${geno}/${pop}
# else
# for ff in {1..4}; do
# mkdir -p /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/ref_data/ldblk_snplist/${split}/${geno}_fold${ff}/${pop}
# done
# fi
# done
# done
# done

module load R/4.2.0-foss-2020b
Rscript /gpfs/gibbs/pi/zhao/xz674/simulation/a.data_prepare/subsample_GWAS/a.geno_data/c.ld_divide.R
