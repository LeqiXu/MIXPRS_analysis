#!/bin/bash
#SBATCH --mem=2G
#SBATCH -p scavenge,day,week,pi_zhao
#SBATCH -t 0:10:00 
#SBATCH -c 1

module load R/4.3.0-foss-2020b

# --- Config ---
BASE_DIR="/home/yd357/pi_paths/pi_zhao/GWAS_Subsample_JointPRS/v2_PUMAS-EN/single-pop/Lassosum"
PARA_DIR="${BASE_DIR}/codes/base"
CORE_R_SCRIPT="${PARA_DIR}/run_lassosum_single.R"

IFS=',' read -r TRAIT ANCESTRY PUMAS ITER CHR S LAMBDA N < <(sed -n "${SLURM_ARRAY_TASK_ID}p" ${PARA_DIR}/job_params.txt)
for var in TRAIT ANCESTRY PUMAS ITER CHR S LAMBDA N; do
    eval "$var=\$(echo \$$var | tr -d '[:space:]')"
done

# --- File Paths ---
if [ "$PUMAS" == "T" ]; then
    GWAS_RAW_FILE="/home/yd357/pi_paths/pi_zhao/GWAS_Subsample_JointPRS/v2_PUMAS-EN/PUMAS-EN/data/subsampling/PUMA-ensemble.subsampling/EN/${TRAIT}_${ANCESTRY}_inter/${TRAIT}_${ANCESTRY}_inter.gwas.omnibus.ite${ITER}.txt"
    FORMATTED_GWAS_DIR="${BASE_DIR}/results/sub_gwas/PUMAS-ite${ITER}"
    FINAL_OUT_DIR="${BASE_DIR}/results/PRSweights/PUMAS-ite${ITER}/${TRAIT}_${ANCESTRY}"
else
    GWAS_RAW_FILE="/gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/data/summary_data/clean/${TRAIT}_${ANCESTRY}_inter_clean.txt"
    FORMATTED_GWAS_DIR="${BASE_DIR}/results/sub_gwas/Regular"
    FINAL_OUT_DIR="${BASE_DIR}/results/PRSweights/Regular/${TRAIT}_${ANCESTRY}"
fi

FORMATTED_GWAS_FILE="${FORMATTED_GWAS_DIR}/${TRAIT}_${ANCESTRY}_chr${CHR}.lassosum_input.txt"
OUT_PREFIX="${FINAL_OUT_DIR}/${TRAIT}_${ANCESTRY}_chr${CHR}"

LDREF_PREFIX="/home/yd357/pi_paths/pi_zhao/GWAS_Subsample_JointPRS/v2_PUMAS-EN/PUMAS-EN/data/LD_ref/A1A2.match/2_make_LD_genotype/GLGC/${ANCESTRY}/qced/1KG.LD500.PRSCS.1kg.A1A2match"
LD_BLOCK_BASE_DIR="/vast/palmer/home.mccleary/yd357/R/x86_64-pc-linux-gnu-library/4.3/lassosum/data"

case "$ANCESTRY" in
    EUR) LD_BLOCK_FILE="${LD_BLOCK_BASE_DIR}/Berisa.EUR.hg19.bed" ;;
    AFR) LD_BLOCK_FILE="${LD_BLOCK_BASE_DIR}/Berisa.AFR.hg19.bed" ;;
    EAS) LD_BLOCK_FILE="${LD_BLOCK_BASE_DIR}/Berisa.ASN.hg19.bed" ;;
    AMR|SAS) LD_BLOCK_FILE="${LD_BLOCK_BASE_DIR}/Berisa.EUR.hg19.bed" ;;
    *) LD_BLOCK_FILE="" ;;
esac

mkdir -p "${FORMATTED_GWAS_DIR}" "${FINAL_OUT_DIR}"

# --- Format GWAS if missing ---
if [ ! -s "${FORMATTED_GWAS_FILE}" ]; then
    awk -v chr_target="${CHR}" 'BEGIN {
        FS="[\t ,]+"; OFS="\t"
    }
    NR==1 {
        for (i=1; i<=NF; i++) col[tolower($i)] = i
        print "CHR", "SNP", "A1", "A2", "BETA", "P", "N"
        next
    }
    $col["chr"] == chr_target {
        print $col["chr"], $col["snp"], $col["a1"], $col["a2"], $col["beta"], $col["p"], $col["n"]
    }' "${GWAS_RAW_FILE}" > "${FORMATTED_GWAS_FILE}"
fi

# --- Run Lassosum ---
Rscript "${CORE_R_SCRIPT}" \
    --gwas_file "${FORMATTED_GWAS_FILE}" \
    --sample_size "${N}" \
    --ref_bfile "${LDREF_PREFIX}" \
    --ld_blocks "${LD_BLOCK_FILE}" \
    --s_param "${S}" \
    --lambda "${LAMBDA}" \
    --out_prefix "${OUT_PREFIX}"
