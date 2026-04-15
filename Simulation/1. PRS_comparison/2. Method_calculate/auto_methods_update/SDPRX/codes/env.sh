#!/bin/bash
#SBATCH --job-name=SDPRX_arreeeeeeeay%j
#SBATCH --output=/home/yd357/pi_paths/pi_zhao/MIXPRS_Revise/Project1_15K_80K/code/4_methods/SDPRX/logs/eneneneneSDPRX_array%j.log
#SBATCH --error=/home/yd357/pi_paths/pi_zhao/MIXPRS_Revise/Project1_15K_80K/code/4_methods/SDPRX/logs/eneeenneneSDPRX_array%j.log
#SBATCH --requeue
#SBATCH --mem=32G
#SBATCH -p scavenge,day,pi_zhao
#SBATCH -t 1:00:00
#SBATCH -c 1

module load miniconda
# Assuming SDPRX is in the same environment as JointPRS
source /vast/palmer/apps/avx2/software/miniconda/24.9.2/etc/profile.d/conda.sh
conda activate /home/yd357/pi_paths/pi_zhao/.conda/envs/JointPRS 

conda install pandas=1.5.3
conda list pandas
