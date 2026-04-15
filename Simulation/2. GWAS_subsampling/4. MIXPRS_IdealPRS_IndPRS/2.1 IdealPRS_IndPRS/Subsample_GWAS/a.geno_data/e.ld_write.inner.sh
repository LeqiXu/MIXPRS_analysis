#!/bin/bash
#SBATCH --partition=scavenge,day
#SBATCH --requeue
#SBATCH --mem=15G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1 --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --job-name=ld_write
#SBATCH --output=/gpfs/gibbs/pi/zhao/xz674/logs/subsample_GWAS/0923ld_write_2/ld_write_%A_%a.txt

params_file="/gpfs/gibbs/pi/zhao/xz674/simulation/a.data_prepare/subsample_GWAS/a.geno_data/d.params.txt"
IFS=' ' read pop split geno ff <<< $(sed -n "${SLURM_ARRAY_TASK_ID}p" ${params_file})

module load miniconda
conda activate py_env

if [[ "${geno}" == "1kg" || "${geno}" == "1kg1" || "${geno}" == "1kg2" ]]; then
    geno_type="1kg"
else
    geno_type="${geno}"
fi

if [[ "${split}" == "discover_validate" ]]; then
mkdir -p /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/ref_data/${split}/${geno}/ldblk_${geno_type}_${pop,,}
else
mkdir -p /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/ref_data/${split}/${geno}_fold${ff}/ldblk_${geno_type}_${pop,,}
fi

python /gpfs/gibbs/pi/zhao/xz674/simulation/a.data_prepare/subsample_GWAS/a.geno_data/e.ld_write.py \
--pop=${pop} \
--split=${split} \
--geno=${geno} \
--fold=${ff}
