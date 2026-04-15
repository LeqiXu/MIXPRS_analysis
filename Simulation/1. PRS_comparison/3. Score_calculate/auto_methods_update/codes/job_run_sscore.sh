#!/bin/bash
#SBATCH --job-name=sscore_array%j
#SBATCH --output=/home/yd357/pi_paths/pi_zhao/MIXPRS_Revise/Project1_15K_80K/code/5_evaluation/1_sscore/logs/25.11.8.2/sscore_array_%A_%a.log
#SBATCH --requeue
#SBATCH --mem=16G
#SBATCH -p scavenge,day,pi_zhao
#SBATCH -t 1:00:00
#SBATCH -c 1

set -e

module load PLINK/2

# --- Helper function to remove duplicate SNPs ---
dedup_beta() {
    local beta_file="$1"
    local dedup_file="${beta_file}.dedup"
    awk '!seen[$1]++' "$beta_file" > "$dedup_file"
    echo "$dedup_file"
}

# --- SLURM Parameters ---
CODE_DIR="/home/yd357/pi_paths/pi_zhao/MIXPRS_Revise/Project1_15K_80K/code/5_evaluation/1_sscore/codes"
IFS=',' read -r METHOD SIM_I P SAMPLE_SIZE POP2 < <(sed -n "${SLURM_ARRAY_TASK_ID}p" ${CODE_DIR}/sscore_params.txt)

# --- Static Parameters ---
RHOG="0.8"
POP1="EUR"

# --- Common Input Dirs ---
GENO_DIR="/gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/data/sim_data/geno_data"
SNPLIST_DIR="/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/snplist_data/simulation"

# --- Input Beta Directories (from 'bind_chromosome' step) ---
BETA_DIR_BASE="/home/yd357/pi_paths/pi_zhao/MIXPRS_Revise/Project1_15K_80K/result/add_benchmark/25.10.4/sim_data/bind_chromosome"
BETA_DIR_JOINTPRS="${BETA_DIR_BASE}/JointPRS"
BETA_DIR_SDPRX="${BETA_DIR_BASE}/SDPRX"
BETA_DIR_XPASS="${BETA_DIR_BASE}/XPASS"
BETA_DIR_MIXPRS_J="${BETA_DIR_BASE}/MIXPRS/MIXPRS-JointPRS"
BETA_DIR_MIXPRS_S="${BETA_DIR_BASE}/MIXPRS/MIXPRS-SDPRX"

# --- Output Score Directory ---
SCORE_DIR_BASE="/home/yd357/pi_paths/pi_zhao/MIXPRS_Revise/Project1_15K_80K/result/add_benchmark/25.10.4/sim_data/prs_scores"

# --- Define Common Input Files for PLINK ---
BFILE_PATH="${GENO_DIR}/${POP2}/test/${POP2}"
SNPLIST_PATH="${SNPLIST_DIR}/${POP2}_sim${SIM_I}_p${P}_rho${RHOG}_${SAMPLE_SIZE}_inter_snplist_real.txt"

if [ ! -f "${BFILE_PATH}.bed" ]; then
    echo "ERROR: Genotype file not found: ${BFILE_PATH}.bed"
    exit 1
fi
if [ ! -f "$SNPLIST_PATH" ]; then
    echo "ERROR: SNP list file not found: $SNPLIST_PATH"
    exit 1
fi

echo "========================================================"
echo "RUNNING: $METHOD"
echo "  Simulation ID : ${SIM_I}"
echo "  P-value       : ${P}"
echo "  Sample Scenario : ${SAMPLE_SIZE}"
echo "  Test Pop (bfile): ${POP2}"
echo "========================================================"

# --- Main Logic ---

if [ "$METHOD" == "JointPRS" ] || [ "$METHOD" == "MIXPRS_J" ]; then
    
    # Set paths based on method
    if [ "$METHOD" == "JointPRS" ]; then
        OUT_DIR="${SCORE_DIR_BASE}/JointPRS"
        BETA_DIR="${BETA_DIR_JOINTPRS}"
    else
        OUT_DIR="${SCORE_DIR_BASE}/MIXPRS/MIXPRS-JointPRS"
        BETA_DIR="${BETA_DIR_MIXPRS_J}"
    fi
    mkdir -p $OUT_DIR
    
    echo "  Generating 5 scores for ${METHOD}..."
    for SCORE_POP in "EUR" "EAS" "AFR" "SAS" "AMR"; do
        BETA_FILE="${BETA_DIR}/sim${SIM_I}_p${P}_rho${RHOG}_${SAMPLE_SIZE}_JointPRS_beta_${SCORE_POP}.txt"
        OUT_NAME="sim${SIM_I}_p${P}_rho${RHOG}_${SAMPLE_SIZE}_JointPRS_beta_${SCORE_POP}_prs_${POP2}"
        
        if [ ! -f "$BETA_FILE" ]; then echo "WARNING: Missing beta file: $BETA_FILE"; continue; fi
        
        BETA_FILE=$(dedup_beta "$BETA_FILE")
        
        plink2 \
          --bfile $BFILE_PATH \
          --extract $SNPLIST_PATH \
          --score $BETA_FILE 1 2 3 header \
          --double-id \
          --threads 1 \
          --out $OUT_DIR/$OUT_NAME
        
        echo "    - Created $OUT_NAME.sscore"
    done

elif [ "$METHOD" == "SDPRX" ] || [ "$METHOD" == "MIXPRS_S" ]; then

    # Set paths based on method
    if [ "$METHOD" == "SDPRX" ]; then
        OUT_DIR="${SCORE_DIR_BASE}/SDPRX"
        BETA_DIR="${BETA_DIR_SDPRX}"
    else
        OUT_DIR="${SCORE_DIR_BASE}/MIXPRS/MIXPRS-SDPRX"
        BETA_DIR="${BETA_DIR_MIXPRS_S}"
    fi
    mkdir -p $OUT_DIR

    echo "  Generating 2 scores for ${METHOD}..."
    for SCORE_POP in "$POP1" "$POP2"; do
        BETA_FILE="${BETA_DIR}/sim${SIM_I}_p${P}_rho${RHOG}_${SAMPLE_SIZE}_${POP1}_${POP2}_beta_${SCORE_POP}.txt"
        OUT_NAME="sim${SIM_I}_p${P}_rho${RHOG}_${SAMPLE_SIZE}_${POP1}_${POP2}_beta_${SCORE_POP}_prs_${POP2}"
        
        if [ ! -f "$BETA_FILE" ]; then echo "WARNING: Missing beta file: $BETA_FILE"; continue; fi

        BETA_FILE=$(dedup_beta "$BETA_FILE")

        plink2 \
          --bfile $BFILE_PATH \
          --extract $SNPLIST_PATH \
          --score $BETA_FILE 1 2 3 header \
          --double-id \
          --threads 1 \
          --out $OUT_DIR/$OUT_NAME
          
        echo "    - Created $OUT_NAME.sscore"
    done

elif [ "$METHOD" == "XPASS" ]; then
    OUT_DIR="${SCORE_DIR_BASE}/XPASS"
    BETA_DIR="${BETA_DIR_XPASS}"
    mkdir -p $OUT_DIR
    
    echo "  Generating 1 score for XPASS..."
    SCORE_POP=$POP2
    BETA_FILE="${BETA_DIR}/sim${SIM_I}_p${P}_rho${RHOG}_${SAMPLE_SIZE}_XPASS_${POP1}_${POP2}_beta_${SCORE_POP}.txt"
    OUT_NAME="sim${SIM_I}_p${P}_rho${RHOG}_${SAMPLE_SIZE}_XPASS_${POP1}_${POP2}_beta_${SCORE_POP}_prs_${POP2}"
    
    if [ ! -f "$BETA_FILE" ]; then echo "ERROR: Missing beta file: $BETA_FILE"; exit 1; fi
    
    BETA_FILE=$(dedup_beta "$BETA_FILE")
    
    plink2 \
      --bfile $BFILE_PATH \
      --extract $SNPLIST_PATH \
      --score $BETA_FILE 1 2 3 header \
      --double-id \
      --threads 1 \
      --out $OUT_DIR/$OUT_NAME

    echo "    - Created $OUT_NAME.sscore"
else
    echo "ERROR: Unknown method specified: $METHOD"
    exit 1
fi

echo "--- Job complete for $METHOD, sim $SIM_I, p $P, sample $SAMPLE_SIZE, pop $POP2 ---"