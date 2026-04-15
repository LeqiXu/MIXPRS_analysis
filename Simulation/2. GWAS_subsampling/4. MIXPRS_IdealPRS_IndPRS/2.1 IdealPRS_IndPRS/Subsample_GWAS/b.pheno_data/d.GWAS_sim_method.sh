#!/bin/bash
#SBATCH --partition=scavenge,day
#SBATCH --requeue
#SBATCH --mem=50G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1 --nodes=1
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --job-name=GWAS_sim_method
#SBATCH --output=/gpfs/gibbs/pi/zhao/xz674/logs/subsample_GWAS/1005GWAS_sim_clean/GWAS_sim_method_EUR_%A.txt

## make file directories
# for pop in EUR EAS AFR SAS AMR; do
#   mkdir -p /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/summary_data/${pop}/discover/clean
#   # mkdir -p /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/summary_data/${pop}/validate/clean
#   # mkdir -p /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/summary_data/${pop}/discover_validate/clean
# done

## generate GWAS for different methods
module load R/4.2.0-foss-2020b

Rscript /gpfs/gibbs/pi/zhao/xz674/simulation/a.data_prepare/subsample_GWAS/b.pheno_data/d.GWAS_sim_method.R EUR

## remove files
# declare -a pops=("AFR" "EUR" "EAS" "SAS" "AMR")
# declare -a splits=("discover" "validate" "discover_validate")
# for pop in "${pops[@]}"; do
# for split in "${splits[@]}"; do
# rm /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/summary_data/${pop}/${split}/${pop}*
# done
# done
