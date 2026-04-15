## Settings
# declare -a pops=("EUR" "AFR" "EAS" "SAS" "AMR")
# declare -a sims=({1..5})  # proportion of causal SNPs
# declare -a h2s=("0.4")  # heritability
# declare -a ps=("0.001" "0.01" "5e-04" "0.1")  # proportion of causal SNPs
# declare -a rhogs=("0.8")  # genetic correlation
# declare -a ffs=({1..4})  # proportion of causal SNPs

# output_file="/gpfs/gibbs/pi/zhao/xz674/simulation/a.data_prepare/subsample_GWAS/b.pheno_data/c.params.txt"
# > $output_file

# for pop in "${pops[@]}"; do
# for i in "${sims[@]}"; do
# for h2 in "${h2s[@]}"; do
# for p in "${ps[@]}"; do
# for rhog in "${rhogs[@]}"; do
# for ff in "${ffs[@]}"; do
# echo "${pop} ${i} ${h2} ${p} ${rhog} ${ff}" >> $output_file
# done
# done
# done
# done
# done
# done

## job array script
params_file="/gpfs/gibbs/pi/zhao/xz674/simulation/a.data_prepare/subsample_GWAS/b.pheno_data/c.params.txt"
job_count=$(wc -l < ${params_file})

sbatch --array=2-80 /gpfs/gibbs/pi/zhao/xz674/simulation/a.data_prepare/subsample_GWAS/b.pheno_data/c.GWAS_sim.inner.sh

# 1-${job_count}