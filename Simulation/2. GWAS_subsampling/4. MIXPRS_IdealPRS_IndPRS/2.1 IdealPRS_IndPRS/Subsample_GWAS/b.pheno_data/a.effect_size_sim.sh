#!/bin/bash
#SBATCH --partition=scavenge,day,week
#SBATCH --requeue
#SBATCH --mem=50G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1 --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --job-name=effect_size_sim
#SBATCH --output=/gpfs/gibbs/pi/zhao/xz674/logs/subsample_GWAS/0922effect_sim/sim_effect_size.txt

module load R/4.2.0-foss-2020b
Rscript /gpfs/gibbs/pi/zhao/xz674/simulation/a.data_prepare/subsample_GWAS/b.pheno_data/a.effect_size_sim.R
