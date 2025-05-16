#!/bin/bash
#SBATCH --mem=8G
#SBATCH -p scavenge,day,week,pi_zhao 
#SBATCH -t 3:00:00
#SBATCH -c 1

# --- Environment Setup ---
module load miniconda
source /vast/palmer/apps/avx2/software/miniconda/24.9.2/etc/profile.d/conda.sh
conda activate JointPRS

# --- Config ---
BASE_DIR="/home/yd357/pi_paths/pi_zhao/GWAS_Subsample_JointPRS/v2_PUMAS-EN/single-pop/PRS-CS"
PARA_DIR="${BASE_DIR}/codes/base"
PYTHON_PATH="/vast/palmer/pi/zhao/yd357/.conda/envs/JointPRS/bin/python3.8"
PRSCS_PATH="/home/yd357/pi_paths/pi_zhao/X-PRS/PRScsx/codes/github/PRScsx.py"

IFS=',' read -r TRAIT ANCESTRY PUMAS ITER CHR PHI N < <(sed -n "${SLURM_ARRAY_TASK_ID}p" ${PARA_DIR}/job_params.txt)
for var in TRAIT ANCESTRY PUMAS ITER CHR PHI N; do
  eval "$var=\$(echo \$$var | tr -d '\"[:space:]')"
done

# --- Path Setup ---
if [ "$PUMAS" == "T" ]; then
  GWAS_RAW="/home/yd357/pi_paths/pi_zhao/GWAS_Subsample_JointPRS/v2_PUMAS-EN/PUMAS-EN/data/subsampling/PUMA-ensemble.subsampling/EN/${TRAIT}_${ANCESTRY}_inter/${TRAIT}_${ANCESTRY}_inter.gwas.omnibus.ite${ITER}.txt"
  SUBSET_GWAS_PATH="${BASE_DIR}/results/sub_gwas/PUMAS-ite${ITER}"
  FINAL_OUT_DIR="${BASE_DIR}/results/PRSweights/PUMAS-ite${ITER}/${TRAIT}_${ANCESTRY}"
else
  GWAS_RAW="/gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/data/summary_data/clean/${TRAIT}_${ANCESTRY}_inter_clean.txt"
  SUBSET_GWAS_PATH="${BASE_DIR}/results/sub_gwas/Regular"
  FINAL_OUT_DIR="${BASE_DIR}/results/PRSweights/Regular/${TRAIT}_${ANCESTRY}"
fi

mkdir -p "${SUBSET_GWAS_PATH}" "${FINAL_OUT_DIR}"

GWAS_OUT="${SUBSET_GWAS_PATH}/${TRAIT}_${ANCESTRY}.prscs.txt"
LDREF="/home/yd357/pi_paths/pi_zhao/GWAS_Subsample_JointPRS/v2_PUMAS-EN/PUMAS-EN/data/LD_ref/A1A2.match/2_make_LD_genotype/GLGC/${ANCESTRY}/qced/1KG.LD500.PRSCS.1kg.A1A2match"
REF_DIR="/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ref_data/PRScsx/1kg"

# --- Format GWAS File for PRS-CS ---
if [ ! -s "${GWAS_OUT}" ] || [ "$(awk 'END{print NR}' "${GWAS_OUT}")" -lt 6 ]; then
  awk 'BEGIN {FS="\t|,| "; OFS="\t"}
       NR==1 {
         for(i=1;i<=NF;i++) {
           if($i=="SNP") snp=i; if($i=="A1") a1=i;
           if($i=="A2") a2=i; if($i=="BETA") beta=i; if($i=="P") p=i
         }
         print "SNP", "A1", "A2", "BETA", "P"
         next
       }
       {print $snp, $a1, $a2, $beta, $p}' "${GWAS_RAW}" > "${GWAS_OUT}"
else
  awk 'NR==1 {print; next} NF==5 {print}' "${GWAS_OUT}" > "${GWAS_OUT}.tmp"
  mv "${GWAS_OUT}.tmp" "${GWAS_OUT}"
fi

# --- Run PRS-CS ---
ARGS=(
  --ref_dir "${REF_DIR}"
  --bim_prefix "${LDREF}"
  --chrom "${CHR}"
  --sst_file "${GWAS_OUT}"
  --n_gwas "${N}"
  --pop "${ANCESTRY}"
  --out_dir "${FINAL_OUT_DIR}"
  --out_name "${TRAIT}"
)
[ "$PHI" != "auto" ] && ARGS+=(--phi "${PHI}")

"${PYTHON_PATH}" "${PRSCS_PATH}" "${ARGS[@]}"
