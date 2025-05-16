#!/bin/bash
#SBATCH --mem=8G
#SBATCH -p scavenge,day,week,pi_zhao
#SBATCH -t 1:00:00
#SBATCH -c 1

module purge
module load R/4.3.0-foss-2020b
conda activate JointPRS

# --- Configuration ---
BASE_DIR="/home/yd357/pi_paths/pi_zhao/GWAS_Subsample_JointPRS/v2_PUMAS-EN/single-pop/LDPred2"
PARA_DIR="${BASE_DIR}/codes/base"
JOB_PARAMS_FILE="${PARA_DIR}/job_params_ldpred2.txt"
LDPRED2_R_SCRIPT="${PARA_DIR}/run_ldpred2_single.R"

IFS=',' read -r TRAIT ANCESTRY PUMAS ITER CHR H2 N MODE P_PARAM H2_RATIO_PARAM SPARSE_PARAM SHRINK_CORR_PARAM \
    < <(sed -n "${SLURM_ARRAY_TASK_ID}p" ${JOB_PARAMS_FILE})

for var in TRAIT ANCESTRY PUMAS ITER CHR H2 N MODE P_PARAM H2_RATIO_PARAM SPARSE_PARAM SHRINK_CORR_PARAM; do
    eval "$var=\$(echo \$$var | xargs)"
done

# --- File Paths ---
if [ "$PUMAS" == "T" ]; then
    RAW_GWAS_FILE="/home/yd357/pi_paths/pi_zhao/GWAS_Subsample_JointPRS/v2_PUMAS-EN/PUMAS-EN/data/subsampling/PUMA-ensemble.subsampling/EN/${TRAIT}_${ANCESTRY}_inter/${TRAIT}_${ANCESTRY}_inter.gwas.omnibus.ite${ITER}.txt"
    FORMATTED_GWAS_DIR="${BASE_DIR}/results/formatted_gwas/PUMAS/ite${ITER}/${TRAIT}_${ANCESTRY}"
    FINAL_OUT_BASE_DIR="${BASE_DIR}/results/PRSweights/PUMAS/ite${ITER}/${TRAIT}_${ANCESTRY}"
else
    RAW_GWAS_FILE="/gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/data/summary_data/clean/${TRAIT}_${ANCESTRY}_inter_clean.txt"
    FORMATTED_GWAS_DIR="${BASE_DIR}/results/formatted_gwas/Regular/${TRAIT}_${ANCESTRY}"
    FINAL_OUT_BASE_DIR="${BASE_DIR}/results/PRSweights/Regular/${TRAIT}_${ANCESTRY}"
fi

FORMATTED_GWAS_FILE="${FORMATTED_GWAS_DIR}/${TRAIT}_${ANCESTRY}_chr${CHR}.formatted.txt"

if [ "$MODE" == "grid" ]; then
    FINAL_OUT_FILE="${FINAL_OUT_BASE_DIR}/${TRAIT}_${ANCESTRY}_chr${CHR}_grid_p${P_PARAM}_h${H2_RATIO_PARAM}_sp${SPARSE_PARAM}.ldpred2.txt"
elif [ "$MODE" == "auto" ]; then
    FINAL_OUT_FILE="${FINAL_OUT_BASE_DIR}/${TRAIT}_${ANCESTRY}_chr${CHR}_auto_sc${SHRINK_CORR_PARAM}.ldpred2.txt"
fi

LDREF_PREFIX="/home/yd357/pi_paths/pi_zhao/GWAS_Subsample_JointPRS/v2_PUMAS-EN/PUMAS-EN/data/LD_ref/A1A2.match/2_make_LD_genotype/GLGC/${ANCESTRY}/qced/1KG.LD500.PRSCS.1kg.A1A2match"

mkdir -p "${FORMATTED_GWAS_DIR}" "${FINAL_OUT_BASE_DIR}"

# --- Format GWAS ---
if [ ! -s "$FORMATTED_GWAS_FILE" ]; then
    awk -v chr_target="$CHR" '
    BEGIN {
        FS = "\t|,| +"; OFS = "\t";
        print "chr\trsid\tpos\ta1\ta0\tmaf\tbeta\tbeta_se\tn_eff"
    }
    NR == 1 {
        for (i = 1; i <= NF; i++) {
            col_name = $i; gsub(/^[[:space:]]+|[[:space:]]+$/, "", col_name);
            header[col_name] = i
        }
        idx_chr = header["CHR"];
        idx_snp = header["SNP"];
        idx_pos = ("POS" in header) ? header["POS"] : (("BP" in header) ? header["BP"] : 0);
        idx_a1 = header["A1"];
        idx_a2 = header["A2"];
        idx_beta = header["BETA"];
        idx_se = header["SE"];
        idx_freq = ("FREQ" in header) ? header["FREQ"] : 0;
        idx_maf = ("MAF" in header) ? header["MAF"] : 0;
        idx_n = ("N" in header) ? header["N"] : 0;
        idx_n_eff = ("N_eff" in header) ? header["N_eff"] : idx_n;
    }
    NR > 1 {
        chr_val = $idx_chr; sub(/^chr/, "", chr_val)
        if (chr_val == chr_target) {
            rsid = $idx_snp;
            pos = $idx_pos;
            a1 = $idx_a1;
            a2 = $idx_a2;
            beta = $idx_beta;
            beta_se = $idx_se;
            n_eff = $idx_n_eff;

            maf_val = -1;
            if (idx_maf > 0 && $idx_maf >= 0 && $idx_maf <= 1) {
                maf_val = $idx_maf
            } else if (idx_freq > 0 && $idx_freq >= 0 && $idx_freq <= 1) {
                freq_a1 = $idx_freq;
                maf_val = (freq_a1 < 0.5) ? freq_a1 : 1 - freq_a1
            }

            if (maf_val >= 0 && maf_val <= 0.5 && beta_se > 0 && n_eff > 0) {
                print chr_val, rsid, pos, a1, a2, maf_val, beta, beta_se, n_eff
            }
        }
    }' "$RAW_GWAS_FILE" > "$FORMATTED_GWAS_FILE"
fi

# --- Run LDPred2 ---
Rscript ${LDPRED2_R_SCRIPT} \
    --gwas="${FORMATTED_GWAS_FILE}" \
    --ld_ref="${LDREF_PREFIX}" \
    --out="${FINAL_OUT_FILE}" \
    --chr="${CHR}" \
    --h2="${H2}" \
    --total_snps="$(($(wc -l < "$RAW_GWAS_FILE") - 1))" \
    --ancestry="${ANCESTRY}" \
    --cores="${SLURM_CPUS_PER_TASK:-1}" \
    --mode="${MODE}" \
    --p="${P_PARAM}" \
    --h2_ratio="${H2_RATIO_PARAM}" \
    --sparse="${SPARSE_PARAM}" \
    --shrink_corr="${SHRINK_CORR_PARAM}"
