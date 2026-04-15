# Step1: Generate parameter list
# input_file="/gpfs/gibbs/pi/zhao/xz674/simulation/b.method_calculate/IndPRS/SDPRX/params2.txt"
# output_file="/gpfs/gibbs/pi/zhao/xz674/simulation/b.method_calculate/IndPRS/popcorn/params2.txt"

# awk '{NF-=2; print}' $input_file | sort -u > $output_file

# Step2: Estimate beta
params_file="/gpfs/gibbs/pi/zhao/xz674/simulation/b.method_calculate/IndPRS/popcorn/params2.txt"
job_count=$(wc -l < ${params_file})

# popcorn
sbatch --array=1-${job_count} /gpfs/gibbs/pi/zhao/xz674/simulation/b.method_calculate/IndPRS/popcorn/popcorn.inner.sh

# 1-${job_count}