#!/bin/bash
#SBATCH --mem=4G
#SBATCH -p scavenge,day,week 
#SBATCH -t 4:00:00
#SBATCH -c 1

module load miniconda
source /vast/palmer/apps/avx2/software/miniconda/24.9.2/etc/profile.d/conda.sh
conda activate base

CODE_DIR="/home/yd357/pi_paths/pi_zhao/GWAS_Subsample_JointPRS/v2_PUMAS-EN/PUMAS-EN/codes/step1_gwas_qc/1.frq_transforming"
INPUT_DIR="/home/yd357/pi_paths/pi_zhao/GWAS_Subsample_JointPRS/data/ref_data/1kg/afreq"
OUTPUT_DIR="/home/yd357/pi_paths/pi_zhao/GWAS_Subsample_JointPRS/Project22_PUMAS-EN_25.2.24/data/frq"

mkdir -p ${OUTPUT_DIR}

# Read parameters
POP=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${CODE_DIR}/job_params.txt)

python ${CODE_DIR}/convert_afreq.py \
    --input_file ${INPUT_DIR}/${POP}.afreq \
    --output_file ${OUTPUT_DIR}/${POP}.frq

echo "Frequency conversion completed for population: ${POP}"