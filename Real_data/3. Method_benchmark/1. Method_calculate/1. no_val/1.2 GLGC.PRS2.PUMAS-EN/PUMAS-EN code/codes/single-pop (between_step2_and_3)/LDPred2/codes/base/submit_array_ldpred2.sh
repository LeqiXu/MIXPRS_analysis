#!/bin/bash
#SBATCH --mem=1G
#SBATCH --partition=scavenge,day,week
#SBATCH -t 0:10:00
#SBATCH -c 1

# --- Configuration ---
BASE_DIR="/home/yd357/pi_paths/pi_zhao/GWAS_Subsample_JointPRS/v2_PUMAS-EN/single-pop/LDPred2"
PARA_DIR="${BASE_DIR}/codes/base"
JOB_PARAMS_FILE="${PARA_DIR}/job_params_ldpred2.txt"
POP_INFO_FILE="/home/yd357/pi_paths/pi_zhao/GWAS_Subsample_JointPRS/v2_PUMAS-EN/PUMAS-EN/data/basic_pop_info/basic_pop_info.csv"
WEIGHTS_BASE_DIR="${BASE_DIR}/results/PRSweights"
MAX_ITERATIONS=4
MAX_JOBS=9990

p_options=("0.001" "0.01" "0.1")
h2_ratios=("0.1" "0.3")
sparse_options=("T" "F")
shrink_corr_auto="0.5"

mkdir -p "${PARA_DIR}"
> "$JOB_PARAMS_FILE"

params_to_run=()

# --- Build Parameter Set ---
{
    read -r _  # Skip header
    while IFS=, read -r trait ancestry N h2; do
        trait=$(echo "$trait" | tr -d '"[:space:]')
        ancestry=$(echo "$ancestry" | tr -d '"[:space:]')
        N=$(echo "$N" | tr -d '"[:space:]')
        h2=$(echo "$h2" | tr -d '"[:space:]')

        [[ ! "$N" =~ ^[0-9]+$ || ! "$h2" =~ ^[0-9]+(\.[0-9]+)?(e-?[0-9]+)?$ ]] && continue

        # --- PUMAS = T ---
        pumas="T"
        raw_gwas_base="/home/yd357/pi_paths/pi_zhao/GWAS_Subsample_JointPRS/v2_PUMAS-EN/PUMAS-EN/data/subsampling/PUMA-ensemble.subsampling/EN/${trait}_${ancestry}_inter"
        out_dir_base="${WEIGHTS_BASE_DIR}/PUMAS"

        if [ -d "$raw_gwas_base" ]; then
            for iter in $(seq 1 $MAX_ITERATIONS); do
                raw_gwas_file="${raw_gwas_base}/${trait}_${ancestry}_inter.gwas.omnibus.ite${iter}.txt"
                [ ! -f "$raw_gwas_file" ] && continue

                for chr in {1..22}; do
                    out_dir="${out_dir_base}/ite${iter}/${trait}_${ancestry}"

                    for p in "${p_options[@]}"; do
                        for h2r in "${h2_ratios[@]}"; do
                            for sp in "${sparse_options[@]}"; do
                                out_file="${out_dir}/${trait}_${ancestry}_chr${chr}_grid_p${p}_h${h2r}_sp${sp}.ldpred2.txt"
                                [ ! -s "$out_file" ] && params_to_run+=("${trait},${ancestry},${pumas},${iter},${chr},${h2},${N},grid,${p},${h2r},${sp},NA")
                            done
                        done
                    done

                    auto_file="${out_dir}/${trait}_${ancestry}_chr${chr}_auto_sc${shrink_corr_auto}.ldpred2.txt"
                    [ ! -s "$auto_file" ] && params_to_run+=("${trait},${ancestry},${pumas},${iter},${chr},${h2},${N},auto,NA,NA,NA,${shrink_corr_auto}")
                done
            done
        fi

        # --- PUMAS = F ---
        pumas="F"
        iter="NA"
        raw_gwas_file="/gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/data/summary_data/clean/${trait}_${ancestry}_inter_clean.txt"
        [ ! -f "$raw_gwas_file" ] && continue

        out_dir_base="${WEIGHTS_BASE_DIR}/Regular"
        for chr in {1..22}; do
            out_dir="${out_dir_base}/${trait}_${ancestry}"

            for p in "${p_options[@]}"; do
                for h2r in "${h2_ratios[@]}"; do
                    for sp in "${sparse_options[@]}"; do
                        out_file="${out_dir}/${trait}_${ancestry}_chr${chr}_grid_p${p}_h${h2r}_sp${sp}.ldpred2.txt"
                        [ ! -s "$out_file" ] && params_to_run+=("${trait},${ancestry},${pumas},${iter},${chr},${h2},${N},grid,${p},${h2r},${sp},NA")
                    done
                done
            done

            auto_file="${out_dir}/${trait}_${ancestry}_chr${chr}_auto_sc${shrink_corr_auto}.ldpred2.txt"
            [ ! -s "$auto_file" ] && params_to_run+=("${trait},${ancestry},${pumas},${iter},${chr},${h2},${N},auto,NA,NA,NA,${shrink_corr_auto}")
        done
    done
} < "$POP_INFO_FILE"

# --- Save Job Parameters ---
printf "%s\n" "${params_to_run[@]}" > "$JOB_PARAMS_FILE"
job_count=${#params_to_run[@]}

# --- Enforce Job Limit ---
if [ "$job_count" -gt "$MAX_JOBS" ]; then
    head -n $MAX_JOBS "$JOB_PARAMS_FILE" > "${JOB_PARAMS_FILE}.tmp"
    mv "${JOB_PARAMS_FILE}.tmp" "$JOB_PARAMS_FILE"
    job_count=$MAX_JOBS
fi

# --- Submit SLURM Job Array ---
n_jobs_to_submit=${1:-$job_count}
[ "$n_jobs_to_submit" -gt 0 ] && sbatch --array=1-${n_jobs_to_submit}%50 "${PARA_DIR}/job_run_ldpred2.sh"
