# Settings
# declare -a pops=("AFR" "EUR" "EAS" "SAS" "AMR")
# declare -a splits=("discover" "discover_validate") #  "validate"
# declare -a genos=("100K") # "1kg" "1kg1" "1kg2"
# declare -a ffs=({1..4})

# output_file="/gpfs/gibbs/pi/zhao/xz674/simulation/a.data_prepare/subsample_GWAS/a.geno_data/d.params.txt"
# > $output_file

# for pop in "${pops[@]}"; do
# for split in "${splits[@]}"; do
# for geno in "${genos[@]}"; do
# if [[ ${geno} == "1kg" || ${geno} == "1kg1" || ${geno} == "1kg2" ]]; then
# geno_type="1kg"
# else
# geno_type="${geno}"
# fi
# if [[ ${split} == "discover_validate" ]]; then
# echo "${pop} ${split} ${geno} 0" >> $output_file
# else
# for ff in "${ffs[@]}"; do
# echo "${pop} ${split} ${geno} ${ff}" >> $output_file
# done
# fi
# done
# done
# done

## job array script
params_file="/gpfs/gibbs/pi/zhao/xz674/simulation/a.data_prepare/subsample_GWAS/a.geno_data/d.params.txt"
job_count=$(wc -l < ${params_file})

sbatch --array=1-${job_count} /gpfs/gibbs/pi/zhao/xz674/simulation/a.data_prepare/subsample_GWAS/a.geno_data/d.ld_calculate.inner.sh
