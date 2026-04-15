# SDPRX
# Step1: Generate parameter list
# declare -a pop2s=("EAS" "AFR" "SAS" "AMR")
# declare -a sample2s=("100K") # ("25K" "90K")
# declare -a sims=({1..5})
# declare -a ps=("0.001" "0.01" "5e-04" "0.1")
# declare -a chrs=({1..22})
# declare -a ffs=({1..4})

# h2=0.4
# rhog=0.8

# pop1="EUR"
# sample1="ukbb"

# params_file="/gpfs/gibbs/pi/zhao/xz674/simulation/b.method_calculate/IndPRS/SDPRX/params.txt"
# > $params_file

# out_dir="/gpfs/gibbs/pi/zhao/xz674/result/sim_result/IndPRS/SDPRX"

# for pop2 in "${pop2s[@]}"; do
# for sample2 in "${sample2s[@]}"; do
# for sim_i in "${sims[@]}"; do
# for p in "${ps[@]}"; do
# for ff in "${ffs[@]}"; do
# for chr in "${chrs[@]}"; do
# subpops=("${pop1}" "${pop2}")
# for subpop in "${subpops[@]}"; do
# # if [[ "${subpop}" != "EUR" ]]; then
# # # EUR_subnonEUR
# # out_name="sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${pop1}_sub${pop2}_${sample1}_${sample2}_SDPRX_real_fold${ff}"
# # else
# # # subEUR_nonEUR
# # out_name="sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_sub${pop1}_${pop2}_${sample1}_${sample2}_SDPRX_real_fold${ff}"
# # fi
# # if [ ! -f "${out_dir}/${out_name}_chr${chr}_1.txt" ] || [ $(stat -c%s "${out_dir}/${out_name}_chr${chr}_1.txt") -lt 102400 ]; then
# echo "${pop1} ${sample1} ${pop2} ${sample2} ${subpop} ${h2} ${sim_i} ${rhog} ${p} ${ff} ${chr}" >> $params_file
# # fi
# done
# done
# done
# done
# done
# done
# done


# Step2: Estimate beta
params_file="/gpfs/gibbs/pi/zhao/xz674/simulation/b.method_calculate/IndPRS/SDPRX/params.txt"
job_count=$(wc -l < ${params_file})

# SDPRX
sbatch --array=1-${job_count} /gpfs/gibbs/pi/zhao/xz674/simulation/b.method_calculate/IndPRS/SDPRX/PRS2.inner.sh

# 1-${job_count}
