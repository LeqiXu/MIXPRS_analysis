#!/bin/bash
#SBATCH --mem=1G
#SBATCH --partition=scavenge,day,week,pi_zhao
#SBATCH -t 0:30:00
#SBATCH -c 1

# --- Configuration ---
BASE_DIR="/home/yd357/pi_paths/pi_zhao/GWAS_Subsample_JointPRS/v2_PUMAS-EN/single-pop/Lassosum"
PARA_DIR="${BASE_DIR}/codes/base"
job_file="${PARA_DIR}/job_params.txt"
POP_INFO="/home/yd357/pi_paths/pi_zhao/GWAS_Subsample_JointPRS/v2_PUMAS-EN/PUMAS-EN/data/basic_pop_info/basic_pop_info.csv"
LOG_DIR="${BASE_DIR}/logs/$(date +%Y.%m.%d)"
S_PARAMS=("0.2" "0.5" "0.9")
LAMBDA_PARAMS=("0.005" "0.01")

mkdir -p "${PARA_DIR}" "${LOG_DIR}"
SLURM_SUBMIT_LOG="${LOG_DIR}/LassoSub_submit-$(date +%Y%m%d)_${SLURM_JOB_ID}.log"
> "$job_file"
exec > >(tee -a ${SLURM_SUBMIT_LOG}) 2>&1 

# --- Helper function for lambda tag ---
get_lambda_tag() {
    case "$1" in
        0.005) echo "lambda_0p005" ;;
        0.01)  echo "lambda_0p01" ;;
        *)     echo "lambda_unexpected_$(echo "$1" | tr -dc '0-9')" ;;
    esac
}

# --- Generate parameter file ---
tail -n +2 "${POP_INFO}" | while IFS=, read -r trait ancestry N h2 || [[ -n "$trait" ]]; do
    trait=$(echo "$trait" | tr -d '"[:space:]')
    ancestry=$(echo "$ancestry" | tr -d '"[:space:]')
    N=$(echo "$N" | tr -d '"[:space:]')

    [[ -z "$trait" || -z "$ancestry" ]] && continue

    # --- PUMAS=T ---
    pumas="T"
    gwas_base="/home/yd357/pi_paths/pi_zhao/GWAS_Subsample_JointPRS/v2_PUMAS-EN/PUMAS-EN/data/subsampling/PUMA-ensemble.subsampling/EN/${trait}_${ancestry}_inter"

    if [ -d "$gwas_base" ]; then
        for iter in {1..4}; do
            gwas_file="${gwas_base}/${trait}_${ancestry}_inter.gwas.omnibus.ite${iter}.txt"
            [ ! -f "$gwas_file" ] && continue

            out_base="${BASE_DIR}/results/PRSweights/PUMAS-ite${iter}/${trait}_${ancestry}"
            for chr in {1..22}; do
                for s in "${S_PARAMS[@]}"; do
                    for lambda in "${LAMBDA_PARAMS[@]}"; do
                        lambda_tag=$(get_lambda_tag "$lambda")
                        out_file="${out_base}/${trait}_${ancestry}_chr${chr}_sX${s}_${lambda_tag}.lassosum.betas.txt"
                        [ ! -s "$out_file" ] && echo "${trait},${ancestry},${pumas},${iter},${chr},${s},${lambda},${N}" >> "$job_file"
                    done
                done
            done
        done
    fi

    # --- PUMAS=F ---
    pumas="F"
    iter=""
    gwas_file="/gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/data/summary_data/clean/${trait}_${ancestry}_inter_clean.txt"
    [ ! -f "$gwas_file" ] && continue

    out_base="${BASE_DIR}/results/PRSweights/Regular/${trait}_${ancestry}"
    for chr in {1..22}; do
        for s in "${S_PARAMS[@]}"; do
            for lambda in "${LAMBDA_PARAMS[@]}"; do
                lambda_tag=$(get_lambda_tag "$lambda")
                out_file="${out_base}/${trait}_${ancestry}_chr${chr}_sX${s}_${lambda_tag}.lassosum.betas.txt"
                [ ! -s "$out_file" ] && echo "${trait},${ancestry},${pumas},${iter},${chr},${s},${lambda},${N}" >> "$job_file"
            done
        done
    done
done

# --- Final Job Count Check ---
job_count=$(wc -l < "$job_file")
max_jobs=9990

if [ "$job_count" -gt "$max_jobs" ]; then
    head -n $max_jobs "$job_file" > "${job_file}.tmp"
    mv "${job_file}.tmp" "$job_file"
    job_count=$max_jobs
fi

n_jobs_to_submit=${1:-$job_count}

# --- Submit SLURM Job ---
if [ "$n_jobs_to_submit" -gt 0 ]; then
    sbatch --output="${BASE_DIR}/logs/LassoRun_array-%A_%a.log" \
           --job-name="LassoRun" \
           --array=1-${n_jobs_to_submit} \
           "${PARA_DIR}/job_run.sh"
fi
