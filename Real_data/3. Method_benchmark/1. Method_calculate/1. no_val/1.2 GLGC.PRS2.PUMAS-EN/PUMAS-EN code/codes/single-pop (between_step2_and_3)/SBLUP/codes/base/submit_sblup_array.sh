#!/bin/bash
#SBATCH --mem=1G
#SBATCH --partition=scavenge,day,week,pi_zhao
#SBATCH -t 0:10:00
#SBATCH -c 1

BASE_DIR="/home/yd357/pi_paths/pi_zhao/GWAS_Subsample_JointPRS/v2_PUMAS-EN/single-pop/SBLUP"
PARA_DIR="${BASE_DIR}/codes/base"
LOG_DIR="${BASE_DIR}/logs"
JOB_PARAM_DIR="${PARA_DIR}/job_params"
job_file="${JOB_PARAM_DIR}/sblup_job_params.txt"

mkdir -p "${PARA_DIR}" "${LOG_DIR}" "${JOB_PARAM_DIR}"
> "$job_file"

POP_INFO="/home/yd357/pi_paths/pi_zhao/GWAS_Subsample_JointPRS/v2_PUMAS-EN/PUMAS-EN/data/basic_pop_info/basic_pop_info.csv"

tail -n +2 "${POP_INFO}" | while IFS=, read -r trait ancestry N h2; do
    trait=$(echo "$trait" | tr -d '"[:space:]')
    ancestry=$(echo "$ancestry" | tr -d '"[:space:]')
    N=$(echo "$N" | tr -d '"[:space:]')
    h2=$(echo "$h2" | tr -d '"[:space:]')

    # --- PUMAS=T ---
    pumas="T"
    gwas_base="/home/yd357/pi_paths/pi_zhao/GWAS_Subsample_JointPRS/v2_PUMAS-EN/PUMAS-EN/data/subsampling/PUMA-ensemble.subsampling/EN/${trait}_${ancestry}_inter"
    for iter in {1..4}; do
        gwas_file="${gwas_base}/${trait}_${ancestry}_inter.gwas.omnibus.ite${iter}.txt"
        [ ! -f "$gwas_file" ] && continue
        for chr in {1..22}; do
            out_file="${BASE_DIR}/results/PRSweights/PUMAS-ite${iter}/${trait}_${ancestry}/chr${chr}.sblup.cojo"
            [ ! -f "$out_file" ] && echo "${trait},${ancestry},${pumas},${iter},${chr},${h2},${N}" >> "$job_file"
        done
    done

    # --- PUMAS=F ---
    pumas="F"
    iter=""
    gwas_file="/gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/data/summary_data/clean/${trait}_${ancestry}_inter_clean.txt"
    [ ! -f "$gwas_file" ] && continue
    for chr in {1..22}; do
        out_file="${BASE_DIR}/results/PRSweights/Regular/${trait}_${ancestry}/chr${chr}.sblup.cojo"
        [ ! -f "$out_file" ] && echo "${trait},${ancestry},${pumas},${iter},${chr},${h2},${N}" >> "$job_file"
    done
done

job_count=$(wc -l < "$job_file")
max_jobs=9990
[ "$job_count" -gt "$max_jobs" ] && head -n "$max_jobs" "$job_file" > "${job_file}.tmp" && mv "${job_file}.tmp" "$job_file" && job_count=$max_jobs

run_script="${PARA_DIR}/run_sblup_job.sh"
[ ! -f "$run_script" ] && exit 0

n_jobs=${1:-$job_count}
[ "$n_jobs" -gt "$job_count" ] && n_jobs=$job_count

[ "$n_jobs" -gt 0 ] && sbatch --array=1-$n_jobs%100 "$run_script"
