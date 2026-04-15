#!/bin/bash
#SBATCH --job-name=Popcorn_array%j
#SBATCH --output=/home/yd357/pi_paths/pi_zhao/MIXPRS_Revise/Project1_15K_80K/code/4_methods/popcorn/logs/25.10.29/Popcorn_array_%A_%a.log
#SBATCH --error=/home/yd357/pi_paths/pi_zhao/MIXPRS_Revise/Project1_15K_80K/code/4_methods/popcorn/logs/25.10.29/Popcorn_array_%A_%a.log
#SBATCH --requeue
#SBATCH --mem=32G
#SBATCH -p scavenge,day,pi_zhao
#SBATCH -t 24:00:00
#SBATCH -c 1

module load miniconda
source /vast/palmer/apps/avx2/software/miniconda/24.9.2/etc/profile.d/conda.sh
conda activate /home/yd357/pi_paths/pi_zhao/MIXPRS_Revise/env/popcorn

CODE_DIR="/home/yd357/pi_paths/pi_zhao/MIXPRS_Revise/Project1_15K_80K/code/4_methods/popcorn/codes"
INPUT_DIR="/home/yd357/pi_paths/pi_zhao/MIXPRS_Revise/Project1_15K_80K/data/sim_data/summary_data/discover_validate/popcorn"
RESULT_DIR="/home/yd357/pi_paths/pi_zhao/MIXPRS_Revise/Project1_15K_80K/result/add_benchmark/25.10.4/sim_data/popcorn"
# Path to Popcorn reference data (cscore files)
POPCORN_REF_DIR="/gpfs/gibbs/pi/zhao/lx94/JointPRS/method/Popcorn/ref"

PARAM_FILE="${CODE_DIR}/popcorn_params.txt"
IFS=',' read -r SIM_I P SAMPLE_SIZE POP2 < <(sed -n "${SLURM_ARRAY_TASK_ID}p" $PARAM_FILE)

RHOG="0.8"
POP1="EUR"
SAMPLE1_LABEL="80K" # EUR (pop1) always uses 80K

SFILE1="${INPUT_DIR}/${POP1}_sim${SIM_I}_p${P}_rho${RHOG}_${SAMPLE1_LABEL}_subcol.txt"
SFILE2="${INPUT_DIR}/${POP2}_sim${SIM_I}_p${P}_rho${RHOG}_${SAMPLE_SIZE}_subcol.txt"
CFILE="${POPCORN_REF_DIR}/${POP1}_${POP2}_all_gen_eff.cscore"
OUT_FILE="${RESULT_DIR}/sim${SIM_I}_p${P}_rho${RHOG}_${POP1}_${POP2}_${SAMPLE1_LABEL}_${SAMPLE_SIZE}_popcorn_corr.txt"

echo "========================================================"
echo "RUNNING Popcorn with parameters:"
echo "  Job ID          : ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
echo "  Simulation ID   : ${SIM_I}"
echo "  P-value         : ${P}"
echo "  Sample Scenario : ${SAMPLE_SIZE} (for ${POP2})"
echo "  Population 1    : ${POP1} (${SAMPLE1_LABEL})"
echo "  Population 2    : ${POP2} (${SAMPLE_SIZE})"
echo "  sfile1 (Pop1)   : ${SFILE1}"
echo "  sfile2 (Pop2)   : ${SFILE2}"
echo "  cfile           : ${CFILE}"
echo "  Output File     : ${OUT_FILE}"
echo "========================================================"

popcorn fit -v 0 \
    --cfile ${CFILE} \
    --gen_effect \
    --sfile1 ${SFILE1} \
    --sfile2 ${SFILE2} \
    ${OUT_FILE}

echo "Popcorn analysis complete for sim${SIM_I}, p=${P}, ${POP1}(${SAMPLE1_LABEL}) vs ${POP2}(${SAMPLE_SIZE})"
echo "Output saved to ${OUT_FILE}"