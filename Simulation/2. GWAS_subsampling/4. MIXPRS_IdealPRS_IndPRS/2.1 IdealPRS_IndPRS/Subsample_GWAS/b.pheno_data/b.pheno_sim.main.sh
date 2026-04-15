## Settings
# declare -a pops=("EUR" "EAS" "AFR" "SAS" "AMR")
# declare -a sims=({1..5})  # proportion of causal SNPs
# declare -a h2s=("0.4")  # heritability
# declare -a ps=("0.001" "0.01" "5e-04" "0.1")  # proportion of causal SNPs
# declare -a rhogs=("0.8")  # genetic correlation

# output_file="/gpfs/gibbs/pi/zhao/xz674/simulation/a.data_prepare/subsample_GWAS/b.pheno_data/b.params.txt"
# > $output_file

# for pop in "${pops[@]}"; do
# for i in "${sims[@]}"; do
# for h2 in "${h2s[@]}"; do
# for p in "${ps[@]}"; do
# for rhog in "${rhogs[@]}"; do
# echo "${pop} ${i} ${h2} ${p} ${rhog}" >> $output_file
# done
# done
# done
# done
# done

## job array script
params_file="/gpfs/gibbs/pi/zhao/xz674/simulation/a.data_prepare/subsample_GWAS/b.pheno_data/b.params.txt"
job_count=$(wc -l < ${params_file})

sbatch --array=4,8,12,16,20 /gpfs/gibbs/pi/zhao/xz674/simulation/a.data_prepare/subsample_GWAS/b.pheno_data/b.pheno_sim.inner.sh


## remove files
# for pop in EUR EAS AFR SAS AMR; do
# for split in discover_validate discover validate test; do
#   rm /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/summary_data/${pop}/${split}/clean/*1kg*
#   rm /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/summary_data/${pop}/${split}/PRScsx/*1kg*
#   rm /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/pheno_data/${pop}/${split}/*1kg*
# done
# done