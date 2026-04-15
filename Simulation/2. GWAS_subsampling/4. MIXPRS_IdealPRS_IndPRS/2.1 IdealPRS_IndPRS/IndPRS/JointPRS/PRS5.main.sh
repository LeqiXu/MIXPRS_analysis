# JointPRS
# Step1: Generate parameter list
# declare -a subpops=("EUR" "EAS" "AFR" "SAS" "AMR") # 
# declare -a sample2s=("100K") # ("25K" "90K")
# declare -a sims=({1..5})
# declare -a ps=("0.001" "0.01" "5e-04" "0.1")
# declare -a chrs=({1..22})
# declare -a ffs=({1..4})

# h2=0.4
# rhog=0.8

# sample1="ukbb"

# params_file="/gpfs/gibbs/pi/zhao/xz674/simulation/b.method_calculate/IndPRS/JointPRS/params.txt"
# > $params_file

# for subpop in "${subpops[@]}"; do
# for sample2 in "${sample2s[@]}"; do
# for sim_i in "${sims[@]}"; do
# for p in "${ps[@]}"; do
# for ff in "${ffs[@]}"; do
# for chr in "${chrs[@]}"; do
# echo "${sample1} ${subpop} ${sample2} ${h2} ${sim_i} ${rhog} ${p} ${ff} ${chr}" >> $params_file
# done
# done
# done
# done
# done
# done

# Step2: Estimate beta
params_file="/gpfs/gibbs/pi/zhao/xz674/simulation/b.method_calculate/IndPRS/JointPRS/params.txt"
job_count=$(wc -l < ${params_file})

# JointPRS
sbatch --array=3301-${job_count} /gpfs/gibbs/pi/zhao/xz674/simulation/b.method_calculate/IndPRS/JointPRS/PRS5.inner.sh

# 1-${job_count}


