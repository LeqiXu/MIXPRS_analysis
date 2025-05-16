#!/bin/bash
#SBATCH --mem=1G
#SBATCH --partition=scavenge,day,week,pi_zhao
#SBATCH -t 0:10:00
#SBATCH -c 1

# --- Configuration ---
BASE_DIR="/home/yd357/pi_paths/pi_zhao/GWAS_Subsample_JointPRS/v2_PUMAS-EN/single-pop/PRS-CS"
PARA_DIR="${BASE_DIR}/codes/base"
POP_INFO="/home/yd357/pi_paths/pi_zhao/GWAS_Subsample_JointPRS/v2_PUMAS-EN/PUMAS-EN/data/basic_pop_info/basic_pop_info.csv"
job_file="${PARA_DIR}/job_params.txt"

mkdir -p "${PARA_DIR}"
> "$job_file"

# --- Generate Parameters ---
tail -n +2 "$POP_INFO" | while IFS=, read -r trait ancestry N h2; do
    trait=$(echo "$trait" | tr -d '"')
    ancestry=$(echo "$ancestry" | tr -d '"')
    N=$(echo "$N" | tr -d '"')

    # PUMAS=T
    gwas_dir="${BASE_DIR}/../../PUMAS-EN/data/subsampling/PUMA-ensemble.subsampling/EN/${trait}_${ancestry}_inter"
    if [ -d "$gwas_dir" ]; then
        for iter in {1..4}; do
            gwas_file="${gwas_dir}/${trait}_${ancestry}_inter.gwas.omnibus.ite${iter}.txt"
            [ ! -f "$gwas_file" ] && continue
            for chr in {1..22}; do
                for phi in 1e-4 1e-6 auto; do
                    phi_fmt="${phi/e/E}"
                    out_dir="${BASE_DIR}/results/PRSweights/PUMAS-ite${iter}/${trait}_${ancestry}"
                    [ "$phi_fmt" == "auto" ] && suffix="phiauto" || suffix="phi${phi_fmt}"
                    out_file="${out_dir}/${trait}_${ancestry}_pst_eff_a1_b0.5_${suffix}_chr${chr}.txt"
                    [ ! -s "$out_file" ] && echo "${trait},${ancestry},T,${iter},${chr},${phi_fmt},${N}" >> "$job_file"
                done
            done
        done
    fi

    # PUMAS=F
    gwas_file="/gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/data/summary_data/clean/${trait}_${ancestry}_inter_clean.txt"
    [ ! -f "$gwas_file" ] && continue
    for chr in {1..22}; do
        for phi in 1e-4 1e-6 auto; do
            phi_fmt="${phi/e/E}"
            out_dir="${BASE_DIR}/results/PRSweights/Regular/${trait}_${ancestry}"
            [ "$phi_fmt" == "auto" ] && suffix="phiauto" || suffix="phi${phi_fmt}"
            out_file="${out_dir}/${trait}_${ancestry}_pst_eff_a1_b0.5_${suffix}_chr${chr}.txt"
            [ ! -s "$out_file" ] && echo "${trait},${ancestry},F,,${chr},${phi_fmt},${N}" >> "$job_file"
        done
    done
done

# --- Submit Jobs ---
job_count=$(wc -l < "$job_file")
max_jobs=9990
[ "$job_count" -gt "$max_jobs" ] && head -n "$max_jobs" "$job_file" > "${job_file}.tmp" && mv "${job_file}.tmp" "$job_file" && job_count=$max_jobs

n_jobs=${1:-$job_count}
[ "$n_jobs" -gt 0 ] && sbatch --array=1-$n_jobs "${PARA_DIR}/job_run.sh"
