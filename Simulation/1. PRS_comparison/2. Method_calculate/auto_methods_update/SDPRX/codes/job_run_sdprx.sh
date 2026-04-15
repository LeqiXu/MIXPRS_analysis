#!/bin/bash
#SBATCH --job-name=SDPRX_array%j
#SBATCH --output=/home/yd357/pi_paths/pi_zhao/MIXPRS_Revise/Project1_15K_80K/code/4_methods/SDPRX/logs/11.2/SDPRX_array%j.log
#SBATCH --error=/home/yd357/pi_paths/pi_zhao/MIXPRS_Revise/Project1_15K_80K/code/4_methods/SDPRX/logs/11.2/SDPRX_array%j.log
#SBATCH --requeue
#SBATCH --mem=32G
#SBATCH -p scavenge,day,pi_zhao
#SBATCH -t 24:00:00
#SBATCH -c 1

module load miniconda
source /vast/palmer/apps/avx2/software/miniconda/24.9.2/etc/profile.d/conda.sh
conda activate /home/yd357/pi_paths/pi_zhao/.conda/envs/JointPRS 

CODE_DIR="/home/yd357/pi_paths/pi_zhao/MIXPRS_Revise/Project1_15K_80K/code/4_methods/SDPRX/codes"
INPUT_DIR="/home/yd357/pi_paths/pi_zhao/MIXPRS_Revise/Project1_15K_80K/data/sim_data/summary_data/discover_validate/SDPRX"
RESULT_DIR="/home/yd357/pi_paths/pi_zhao/MIXPRS_Revise/Project1_15K_80K/result/add_benchmark/25.10.4/sim_data/SDPRX"
POPCORN_DIR="/home/yd357/pi_paths/pi_zhao/MIXPRS_Revise/Project1_15K_80K/result/add_benchmark/25.10.4/sim_data/popcorn"

SDPRX_PATH="/gpfs/gibbs/pi/zhao/lx94/JointPRS/method/SDPRX/SDPRX.py"
REF_DIR_BASE="/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ref_data/SDPRX"
BIM_PREFIX="/gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/data/sim_data/geno_data/All/All_test"

IFS=',' read -r SIM_I P SAMPLE_SIZE_2 POP2 CHR < <(sed -n "${SLURM_ARRAY_TASK_ID}p" ${CODE_DIR}/sdprx_params.txt)
RHOG="0.8"
POP1="EUR"
SAMPLE_SIZE_1="80K"

N1="311600"
if [[ "${SAMPLE_SIZE_2}" == "15K" ]]; then
    N2="15000"
elif [[ "${SAMPLE_SIZE_2}" == "80K" ]]; then
    N2="80000"
else
    echo "Error: Unknown sample size specified: ${SAMPLE_SIZE_2}"
    exit 1
fi

# --- Get Genetic Correlation (rho_est) ---
POPCORN_FILE="${POPCORN_DIR}/sim${SIM_I}_p${P}_rho${RHOG}_${POP1}_${POP2}_${SAMPLE_SIZE_1}_${SAMPLE_SIZE_2}_popcorn_corr.txt"
if [ ! -f "$POPCORN_FILE" ]; then
    echo "Error: Popcorn file not found: $POPCORN_FILE"
    exit 1
fi
RHO_EST=$(grep '^pge' "${POPCORN_FILE}" | awk '{printf "%.2f", $2}')
if [ -z "$RHO_EST" ]; then
    echo "Error: Could not extract rho_est from $POPCORN_FILE"
    exit 1
fi

# --- Define I/O files ---
LD_DIR="${REF_DIR_BASE}/EUR_${POP2}"
SS1_FILE="${INPUT_DIR}/${POP1}_sim${SIM_I}_p${P}_rho${RHOG}_${SAMPLE_SIZE_1}_subcol.txt"
SS2_FILE="${INPUT_DIR}/${POP2}_sim${SIM_I}_p${P}_rho${RHOG}_${SAMPLE_SIZE_2}_subcol.txt"

# format: sim..._p..._sample2_pop1_pop2_..._chr...
OUT_NAME="sim${SIM_I}_p${P}_rho${RHOG}_${SAMPLE_SIZE_2}_${POP1}_${POP2}_SDPRX_chr${CHR}"

if [ ! -f "$SS1_FILE" ]; then echo "Error: ss1 file not found: $SS1_FILE"; exit 1; fi
if [ ! -f "$SS2_FILE" ]; then echo "Error: ss2 file not found: $SS2_FILE"; exit 1; fi
if [ ! -d "$LD_DIR" ]; then echo "Error: LD reference directory not found: $LD_DIR"; exit 1; fi

mkdir -p "${RESULT_DIR}"

echo "========================================================"
echo "RUNNING SDPRX with parameters:"
echo "  Simulation ID : ${SIM_I}"
echo "  P-value       : ${P}"
echo "  Sample Scenario : ${SAMPLE_SIZE_2} (for ${POP2})"
echo "  Population 1  : ${POP1} (N=${N1})"
echo "  Population 2  : ${POP2} (N=${N2})"
echo "  Chromosome    : ${CHR}"
echo "  Est. Rho      : ${RHO_EST}"
echo "  LD Directory  : ${LD_DIR}"
echo "  SS1 File      : ${SS1_FILE}"
echo "  SS2 File      : ${SS2_FILE}"
echo "  Out Name      : ${OUT_NAME}"
echo "========================================================"

python $SDPRX_PATH \
    --load_ld $LD_DIR \
    --valid ${BIM_PREFIX}.bim \
    --ss1 $SS1_FILE \
    --ss2 $SS2_FILE \
    --N1 $N1 \
    --N2 $N2 \
    --mcmc_samples 2000 \
    --burn 1000 \
    --force_shared True \
    --chr ${CHR} \
    --rho ${RHO_EST} \
    --out ${RESULT_DIR}/${OUT_NAME}

echo "SDPRX analysis complete for sim${SIM_I}, p=${P}, sample_size=${SAMPLE_SIZE_2}, pop2=${POP2}, chr=${CHR}"