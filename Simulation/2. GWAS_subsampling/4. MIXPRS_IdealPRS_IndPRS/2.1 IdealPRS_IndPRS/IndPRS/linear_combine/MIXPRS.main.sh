## Generate params
# Step1: Generate parameter list
# declare -a subpops=("EUR" "EAS" "AFR" "SAS" "AMR") # 
# declare -a sample2s=("100K") # ("25K" "90K")
# declare -a sims=({1..5})
# declare -a ps=("0.001" "0.01" "5e-04" "0.1")
# declare -a models=("linear" "lasso" "ridge" "elasticnet" "nnls")

# h2=0.4
# rhog=0.8

# type="full_snplist"
# type_original="all"

# sample1="UKB"

# params_file="/gpfs/gibbs/pi/zhao/xz674/simulation/b.method_calculate/IndPRS/linear_combine/MIXPRS.params.txt"
# > $params_file

# for subpop in "${subpops[@]}"; do
# for sample2 in "${sample2s[@]}"; do
# for sim_i in "${sims[@]}"; do
# for p in "${ps[@]}"; do
# for model in "${models[@]}"; do
# echo "${subpop} ${sample1} ${sample2} ${h2} ${sim_i} ${rhog} ${p} ${type} ${type_original} ${model}" >> $params_file
# done
# done
# done
# done
# done

# Step2: Obtain final MIXPRS weight
params_file="/gpfs/gibbs/pi/zhao/xz674/simulation/b.method_calculate/IndPRS/linear_combine/MIXPRS.params.txt"
job_count=$(wc -l < ${params_file})

sbatch --array=2-${job_count} /gpfs/gibbs/pi/zhao/xz674/simulation/b.method_calculate/IndPRS/linear_combine/MIXPRS.inner.sh

# 1-${job_count}