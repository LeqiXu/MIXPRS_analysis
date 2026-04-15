#!/bin/bash
#SBATCH --partition=scavenge,day
#SBATCH --requeue
#SBATCH --mem=10G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1 --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --job-name=Final_weight
#SBATCH --output=/gpfs/gibbs/pi/zhao/xz674/logs/IndPRS/1008Final_weight_JointPRS/Final_weight_JointPRS_calc_%A_%a.txt

## Intergrate final weight
module load R/4.2.0-foss-2020b
Rscript /gpfs/gibbs/pi/zhao/xz674/simulation/b.method_calculate/IndPRS/JointPRS/final_weight.R

## Clean the previous result
rm /gpfs/gibbs/pi/zhao/xz674/result/sim_result/IndPRS/JointPRS/*phiauto_chr*.txt

