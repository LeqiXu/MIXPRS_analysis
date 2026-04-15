# Settings
# Use d.params.txt

## job array script
params_file="/gpfs/gibbs/pi/zhao/xz674/simulation/a.data_prepare/subsample_GWAS/a.geno_data/d.params.txt"
job_count=$(wc -l < ${params_file})

sbatch --array=1-${job_count} /gpfs/gibbs/pi/zhao/xz674/simulation/a.data_prepare/subsample_GWAS/a.geno_data/e.ld_write.inner.sh
