#!/bin/bash
#SBATCH --partition=scavenge,day
#SBATCH --requeue
#SBATCH --mem=10G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1 --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --job-name=Final_weight
#SBATCH --output=/gpfs/gibbs/pi/zhao/xz674/logs/IndPRS/1010Final_weight_SDPRX/Final_weight_SDPRX_calc_%A_%a.txt

## Intergrate final weight
module load R/4.2.0-foss-2020b
Rscript /gpfs/gibbs/pi/zhao/xz674/simulation/b.method_calculate/IndPRS/SDPRX/final_weight.R

## Clean the previous result
# rm /gpfs/gibbs/pi/zhao/xz674/result/sim_result/IndPRS/SDPRX/*chr*.txt

