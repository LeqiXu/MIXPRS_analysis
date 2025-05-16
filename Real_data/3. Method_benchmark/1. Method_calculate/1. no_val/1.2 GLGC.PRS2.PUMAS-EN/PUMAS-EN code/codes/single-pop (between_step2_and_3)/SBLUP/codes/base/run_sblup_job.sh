#!/bin/bash
#SBATCH --mem=8G
#SBATCH -p scavenge,day,week,pi_zhao
#SBATCH -t 1:00:00
#SBATCH -c 1

# --- Environment Setup ---
module purge
module load GCTA/1.94.1-gfbf-2022b
module load R/4.3.0-foss-2020b

R_EXEC=$(which Rscript)
GCTA_EXEC=$(which gcta64)

# --- Config ---
BASE_DIR="/home/yd357/pi_paths/pi_zhao/GWAS_Subsample_JointPRS/v2_PUMAS-EN/single-pop/SBLUP"
PARA_DIR="${BASE_DIR}/codes/base"
JOB_PARAM_DIR="${PARA_DIR}/job_params"
job_file="${JOB_PARAM_DIR}/sblup_job_params.txt"
FORMAT_R_SCRIPT="${PARA_DIR}/format_gwas_for_sblup.R"

IFS=',' read -r TRAIT ANCESTRY PUMAS ITER CHR H2 N_OVERALL < <(sed -n "${SLURM_ARRAY_TASK_ID}p" "$job_file")

# --- Path Setup ---
TEMP_DIR="${BASE_DIR}/temp/${TRAIT}_${ANCESTRY}_${PUMAS:-Regular}${ITER:+_ite${ITER}}_chr${CHR}_task${SLURM_ARRAY_TASK_ID}"
mkdir -p "$TEMP_DIR"

if [ "$PUMAS" == "T" ]; then
    GWAS_RAW_FILE="/home/yd357/pi_paths/pi_zhao/GWAS_Subsample_JointPRS/v2_PUMAS-EN/PUMAS-EN/data/subsampling/PUMA-ensemble.subsampling/EN/${TRAIT}_${ANCESTRY}_inter/${TRAIT}_${ANCESTRY}_inter.gwas.omnibus.ite${ITER}.txt"
    FINAL_OUT_DIR="${BASE_DIR}/results/PRSweights/PUMAS-ite${ITER}/${TRAIT}_${ANCESTRY}"
    GWAS_MA_PATH="${BASE_DIR}/results/formatted_gwas/PUMAS-ite${ITER}/${TRAIT}_${ANCESTRY}"
else
    GWAS_RAW_FILE="/gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/data/summary_data/clean/${TRAIT}_${ANCESTRY}_inter_clean.txt"
    FINAL_OUT_DIR="${BASE_DIR}/results/PRSweights/Regular/${TRAIT}_${ANCESTRY}"
    GWAS_MA_PATH="${BASE_DIR}/results/formatted_gwas/Regular/${TRAIT}_${ANCESTRY}"
fi

mkdir -p "$FINAL_OUT_DIR" "$GWAS_MA_PATH"

GWAS_MA_FILE="${GWAS_MA_PATH}/${TRAIT}_${ANCESTRY}_chr${CHR}.ma"
SNP_COUNT_M_FILE="${TEMP_DIR}/snp_count_m.txt"

# --- Format GWAS with R ---
"$R_EXEC" "$FORMAT_R_SCRIPT" "$GWAS_RAW_FILE" "$GWAS_MA_FILE" "$CHR" "$SNP_COUNT_M_FILE"

# --- Compute --cojo-sblup ---
m=$(cat "$SNP_COUNT_M_FILE")
cojo_sblup=$(printf "%.0f" "$(echo "$m * (1 / $H2 - 1)" | bc -l)")

# --- Run GCTA SBLUP ---
LDREF_PREFIX="/home/yd357/pi_paths/pi_zhao/GWAS_Subsample_JointPRS/v2_PUMAS-EN/PUMAS-EN/data/LD_ref/A1A2.match/2_make_LD_genotype/GLGC/${ANCESTRY}/qced/1KG.LD500.PRSCS.1kg.A1A2match"
GCTA_OUT_PREFIX="${FINAL_OUT_DIR}/chr${CHR}"
LD_WINDOW=200
THREAD_NUM=${SLURM_CPUS_PER_TASK:-1}

"$GCTA_EXEC" \
    --bfile "$LDREF_PREFIX" \
    --chr "$CHR" \
    --cojo-file "$GWAS_MA_FILE" \
    --cojo-sblup "$cojo_sblup" \
    --cojo-wind "$LD_WINDOW" \
    --thread-num "$THREAD_NUM" \
    --out "$GCTA_OUT_PREFIX"

# --- Clean Up (optional) ---
# rm -rf "$TEMP_DIR"
