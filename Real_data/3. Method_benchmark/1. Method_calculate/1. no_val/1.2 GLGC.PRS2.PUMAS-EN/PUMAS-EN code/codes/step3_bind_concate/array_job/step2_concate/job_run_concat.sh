#!/bin/bash
#SBATCH --mem=4G
#SBATCH --partition=scavenge,day,week,pi_zhao
#SBATCH -t 3:00:00
#SBATCH -c 1

# --- Configuration ---
BASE_DIR="/home/yd357/pi_paths/pi_zhao/GWAS_Subsample_JointPRS/v3/concatenate"
CODE_DIR="${BASE_DIR}/codes/array_job"
JOB_PARAMS_FILE="${CODE_DIR}/job_params.txt"
CONCAT_BASE_OUTPUT_DIR="${BASE_DIR}/results"
CONCAT_FINAL_DIR="${BASE_DIR}/results/concatenation"
CHROMOSOMES=( $(seq 1 22) )

IFS=',' read -r TRAIT ANCESTRY CASE < <(sed -n "${SLURM_ARRAY_TASK_ID}p" "${JOB_PARAMS_FILE}")
[ -z "${TRAIT}" ] || [ -z "${ANCESTRY}" ] || [ -z "${CASE}" ] && exit 0

trait_anc_dir_name="${TRAIT}_${ANCESTRY}"
trait_anc_file_prefix="${TRAIT}_${ANCESTRY}"

if [ "${CASE}" == "Regular" ]; then
    PUMAS_ITER=""
    case_tag="_Regular"
    out_case_dir_final_output="Regular"
else
    PUMAS_ITER="${CASE#PUMASite}"
    case_tag="_${CASE}"
    out_case_dir_final_output="PUMAS-ite${PUMAS_ITER}"
fi

final_concat_dir="${CONCAT_FINAL_DIR}/${out_case_dir_final_output}/${trait_anc_dir_name}"
mkdir -p "${final_concat_dir}"
final_target_output_file="${final_concat_dir}/${trait_anc_file_prefix}${case_tag}_ALL_METHODS.txt"

TEMP_DIR=$(mktemp -d -p "${SLURM_TMPDIR:-/tmp}" --suffix="_concat_${SLURM_ARRAY_TASK_ID}") 
all_snps_temp_file="${TEMP_DIR}/all_snps.txt"
unique_snps_temp_file="${TEMP_DIR}/unique_snps.txt"
snp_a1_map_temp_file="${TEMP_DIR}/snp_a1_map.txt"
param_files_temp_file="${TEMP_DIR}/param_files.txt"
: > "${all_snps_temp_file}"
: > "${param_files_temp_file}"

# --- Collect SNPs across methods ---
for method in "Lassosum" "LDPred2" "PRS-CS" "SBLUP"; do
    input_dir="${CONCAT_BASE_OUTPUT_DIR}/${method}/${trait_anc_dir_name}"
    [ ! -d "${input_dir}" ] && continue

    pattern="${trait_anc_file_prefix}${case_tag}_ALLCHR_*"
    [ "${method}" == "PRS-CS" ] && pattern="${trait_anc_file_prefix}${case_tag}_*_ALLCHR.txt"

    find "${input_dir}" -maxdepth 1 -type f -name "${pattern}" -print0 | while IFS= read -r -d $'\0' file; do
        case "${method}" in
            "PRS-CS")  awk_cmd='{print $2}' ;;
            *)         awk_cmd='{print $1}' ;;
        esac
        tail -n +2 "$file" | awk "${awk_cmd}" >> "${all_snps_temp_file}"
    done
done

[ ! -s "${all_snps_temp_file}" ] && { echo -e "SNP\tA1" > "${final_target_output_file}"; rm -rf "${TEMP_DIR}"; exit 0; }

LC_ALL=C sort -u "${all_snps_temp_file}" > "${unique_snps_temp_file}"
[ ! -s "${unique_snps_temp_file}" ] && { echo -e "SNP\tA1" > "${final_target_output_file}"; rm -rf "${TEMP_DIR}"; exit 0; }

# --- Select a reliable file for A1 column ---
for method in "Lassosum" "LDPred2" "PRS-CS" "SBLUP"; do
    input_dir="${CONCAT_BASE_OUTPUT_DIR}/${method}/${trait_anc_dir_name}"
    [ ! -d "${input_dir}" ] && continue
    pattern="${trait_anc_file_prefix}${case_tag}_ALLCHR_*"
    [ "${method}" == "PRS-CS" ] && pattern="${trait_anc_file_prefix}${case_tag}_*_ALLCHR.txt"
    file=$(find "${input_dir}" -maxdepth 1 -type f -name "${pattern}" | head -n 1)
    [ -f "$file" ] || continue

    case "${method}" in
        "Lassosum"|"LDPred2") awk_cmd='{print $1 "\t" $2}' ;;
        "PRS-CS")             awk_cmd='{print $2 "\t" $4}' ;;
        "SBLUP")              awk_cmd='{print $1 "\t" $2}' ;;
    esac
    tail -n +2 "$file" | awk "${awk_cmd}" > "${snp_a1_map_temp_file}"
    [ -s "${snp_a1_map_temp_file}" ] && break
done

# --- Prepare param files for R script ---
for method in "Lassosum" "LDPred2" "PRS-CS" "SBLUP"; do
    input_dir="${CONCAT_BASE_OUTPUT_DIR}/${method}/${trait_anc_dir_name}"
    [ ! -d "${input_dir}" ] && continue

    pattern="${trait_anc_file_prefix}${case_tag}_ALLCHR_*"
    [ "${method}" == "PRS-CS" ] && pattern="${trait_anc_file_prefix}${case_tag}_*_ALLCHR.txt"

    find "${input_dir}" -maxdepth 1 -type f -name "${pattern}" -print0 | while IFS= read -r -d $'\0' file; do
        filename=$(basename "${file}")
        params=""
        if [ "${method}" == "PRS-CS" ]; then
            params="${filename#"${trait_anc_file_prefix}${case_tag}_"}"
            params="${params%"_ALLCHR.txt"}"
        elif [ "${method}" == "SBLUP" ]; then
            params="${filename#"${trait_anc_file_prefix}${case_tag}_ALLCHR_"}"
            params="${params%".sblup.cojo"}"
            [ -z "$params" ] && params="beta"
        else
            suffix=".lassosum.betas.txt"
            [ "${method}" == "LDPred2" ] && suffix=".ldpred2.txt"
            params="${filename#"${trait_anc_file_prefix}${case_tag}_ALLCHR_"}"
            params="${params%$suffix}"
        fi
        echo -e "${method}\t${method}_${params}\t${file}" >> "${param_files_temp_file}"
    done
done

# --- Run R script ---
cd "${TEMP_DIR}"
module load R/4.3.0-foss-2020b
Rscript_output=$(Rscript -e '
snps <- read.table("unique_snps.txt", header=FALSE, stringsAsFactors=FALSE); colnames(snps) <- "SNP"
a1 <- tryCatch(read.table("snp_a1_map.txt", header=FALSE), error=function(e) NULL)
if (!is.null(a1) && nrow(a1) > 0) { colnames(a1) <- c("SNP", "A1"); res <- merge(snps, a1, by="SNP", all.x=TRUE) } else { res <- data.frame(SNP=snps$SNP, A1=NA_character_) }
params <- tryCatch(read.table("param_files.txt", header=FALSE, sep="\t", stringsAsFactors=FALSE), error=function(e) NULL)
if (!is.null(params)) { colnames(params) <- c("method", "column_name", "filename")
for (i in 1:nrow(params)) {
  d <- tryCatch({
    if (!file.exists(params$filename[i])) return(NULL)
    f <- read.table(params$filename[i], header=TRUE, stringsAsFactors=FALSE)
    cname <- params$column_name[i]
    if (params$method[i] == "Lassosum" && all(c("SNP","BETA") %in% names(f))) data.frame(SNP=f$SNP, val=f$BETA) else
    if (params$method[i] == "LDPred2" && all(c("rsid","beta") %in% names(f))) data.frame(SNP=f$rsid, val=f$beta) else
    if (params$method[i] == "PRS-CS" && all(c("SNP","BETA") %in% names(f))) data.frame(SNP=f$SNP, val=f$BETA) else
    if (params$method[i] == "SBLUP" && "SNP" %in% names(f)) {
      bcol <- if ("BETA_SBLUP" %in% names(f)) "BETA_SBLUP" else if ("b_SBLUP" %in% names(f)) "b_SBLUP" else names(f)[4]
      data.frame(SNP=f$SNP, val=f[[bcol]])
    } else NULL }, error=function(e) NULL)
  if (!is.null(d)) { names(d) <- c("SNP", params$column_name[i]); res <- merge(res, d, by="SNP", all.x=TRUE) } else {
    res[[params$column_name[i]]] <- NA_real_ }
}}
cols <- setdiff(names(res), c("SNP", "A1"))
for (col in cols) res[[col]][is.na(res[[col]])] <- 0
res <- res[, c("SNP", "A1", sort(cols))]
write.table(res, "final_result.txt", row.names=FALSE, sep="\t", quote=FALSE)
')

# --- Final Output ---
[ -f "${TEMP_DIR}/final_result.txt" ] && cp "${TEMP_DIR}/final_result.txt" "${final_target_output_file}"
rm -rf "${TEMP_DIR}"
