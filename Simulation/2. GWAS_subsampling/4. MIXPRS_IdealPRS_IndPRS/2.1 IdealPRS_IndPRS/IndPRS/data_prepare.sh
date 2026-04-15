#!/bin/bash
#SBATCH --partition=scavenge,day
#SBATCH --requeue
#SBATCH --mem=10G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1 --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --job-name=write_method_data
#SBATCH --output=/gpfs/gibbs/pi/zhao/xz674/logs/IndPRS/1023linear_combine1/write_method_data_%A.txt

module load R/4.2.0-foss-2020b

# Rscript /gpfs/gibbs/pi/zhao/xz674/simulation/b.method_calculate/IndPRS/SDPRX/data_prepare.R

# Rscript /gpfs/gibbs/pi/zhao/xz674/simulation/b.method_calculate/IndPRS/JointPRS/data_prepare.R

# Rscript /gpfs/gibbs/pi/zhao/xz674/simulation/b.method_calculate/IndPRS/popcorn/data_prepare.R

Rscript /gpfs/gibbs/pi/zhao/xz674/simulation/b.method_calculate/IndPRS/linear_combine/PRS.train.R