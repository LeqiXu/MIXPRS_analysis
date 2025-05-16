#!/bin/bash
#SBATCH --mem=4G
#SBATCH --partition=scavenge,day,week,pi_zhao
#SBATCH -t 1:00:00
#SBATCH -c 1

module load R/4.3.0-foss-2020b

# --- Configuration ---
BASE_PROJECT_DIR="/home/yd357/pi_paths/pi_zhao/GWAS_Subsample_JointPRS/v3/concatenate"
CODE_DIR_STEP3="${BASE_PROJECT_DIR}/codes/array_job/step3_gather_trait"
JOB_PARAMS_FILE="${CODE_DIR_STEP3}/job_params_step3_singlefile.txt"
CONCAT_RESULTS_BASE_DIR="${BASE_PROJECT_DIR}/results/concatenation"
SUBSAMPLING_BASE_DIR="/home/yd357/pi_paths/pi_zhao/GWAS_Subsample_JointPRS/v2_PUMAS-EN/PUMAS-EN/data/subsampling/PUMA-ensemble.subsampling"
SUBSAMPLING_ENSEMBLE_METHOD_FOLDER="EN"
FINAL_OUTPUT_BASE_DIR="${BASE_PROJECT_DIR}/results/concatenation_trait_single_file"
POP_INFO_FILE="/home/yd357/pi_paths/pi_zhao/GWAS_Subsample_JointPRS/v2_PUMAS-EN/PUMAS-EN/data/basic_pop_info/basic_pop_info.csv"
ORIGINAL_METHODS_LIST=("Lassosum" "LDPred2" "PRS-CS" "SBLUP")

# --- Read parameters from job file ---
IFS=',' read -r TARGET_TRAIT TARGET_POPULATION_SHELL CASE_TYPE < <(sed -n "${SLURM_ARRAY_TASK_ID}p" "${JOB_PARAMS_FILE}")
[ -z "${TARGET_TRAIT}" ] || [ -z "${TARGET_POPULATION_SHELL}" ] || [ -z "${CASE_TYPE}" ] && exit 1

# --- Case and Path Setup ---
CASE_TAG="${CASE_TYPE}"
if [[ "${CASE_TYPE}" == "Regular" ]]; then
    CASE_DIR_NAME_CONCAT="Regular"
    CASE_TAG_CONCAT="_Regular"
else
    ITER_NUM="${CASE_TYPE#PUMASite}"
    CASE_DIR_NAME_CONCAT="PUMAS-ite${ITER_NUM}"
    CASE_TAG_CONCAT="_${CASE_TYPE}"
fi

PUMAS_TRAIT_NAME_BASE="${TARGET_TRAIT}_${TARGET_POPULATION_SHELL}_${CASE_TAG}"
FINAL_OUTPUT_FILENAME="${PUMAS_TRAIT_NAME_BASE}.all_sources_all_methods.txt"
FINAL_OUTPUT_FILE_PATH="${FINAL_OUTPUT_BASE_DIR}/${FINAL_OUTPUT_FILENAME}"
REF_XTY_FILE="${SUBSAMPLING_BASE_DIR}/${SUBSAMPLING_ENSEMBLE_METHOD_FOLDER}/${TARGET_TRAIT}_${TARGET_POPULATION_SHELL}_inter/${TARGET_TRAIT}_${TARGET_POPULATION_SHELL}_inter.xty.omnibus.ite1.txt"

[ ! -f "${REF_XTY_FILE}" ] && exit 1

mkdir -p "${FINAL_OUTPUT_BASE_DIR}"

# --- Temporary workspace ---
TEMP_DIR=$(mktemp -d -p "${SLURM_TMPDIR:-/tmp}" --suffix="_gatherSF_${SLURM_ARRAY_TASK_ID}")
SOURCE_POP_LIST_FILE="${TEMP_DIR}/source_pop_list.txt"
awk -F, -v tt="${TARGET_TRAIT}" '$1 ~ "\"" tt "\"" {gsub(/"/, "", $2); print $2}' "${POP_INFO_FILE}" | sort -u > "${SOURCE_POP_LIST_FILE}"

# --- Generate embedded R script ---
R_SCRIPT_PATH="${TEMP_DIR}/process_trait_data_singlefile.R"
cat > "${R_SCRIPT_PATH}" << 'EOF_RSCRIPT'
# --- Gather Trait Data into Single File ---
suppressMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)

target_trait <- args[1]
target_pop <- args[2]
case_dir <- args[3]
case_tag <- args[4]
ref_xty <- args[5]
concat_base <- args[6]
final_out <- args[7]
source_pop_list_file <- args[8]
methods <- strsplit(args[9], ",")[[1]]

ref <- fread(ref_xty, select = c("CHR", "SNP", "A1", "A2"), col.names = c("CHR", "SNP", "XTY_A1", "XTY_A2"))
master <- copy(ref)

if (file.exists(source_pop_list_file)) {
  pops <- fread(source_pop_list_file, header = FALSE)$V1
  for (pop in pops) {
    concat_file <- file.path(concat_base, case_dir, paste0(target_trait, "_", pop),
                    paste0(target_trait, "_", pop, case_tag, "_ALL_METHODS.txt"))
    if (!file.exists(concat_file)) next
    dt <- fread(concat_file)
    if (!"SNP" %in% names(dt) || !"A1" %in% names(dt)) next
    setnames(dt, "A1", "Source_A1", skip_absent=TRUE)
    setkey(dt, SNP)

    for (m in methods) {
      cols <- grep(paste0("^", m, "_"), names(dt), value = TRUE)
      for (c in cols) {
        new_col <- paste(pop, c, sep = "_")
        beta_dt <- dt[!is.na(get(c)), .(SNP, Source_A1, BETA = get(c))]
        if (nrow(beta_dt) == 0) {
          master[, (new_col) := 0]
          next
        }
        setkey(beta_dt, SNP)
        merged <- merge(master[, .(SNP, XTY_A1)], beta_dt, by = "SNP", all.x = TRUE)
        merged <- merged[match(master$SNP, merged$SNP)]
        merged[!is.na(Source_A1) & !is.na(XTY_A1) & toupper(Source_A1) != toupper(XTY_A1), BETA := -BETA]
        merged[is.na(BETA), BETA := 0]
        master[, (new_col) := merged$BETA]
      }
    }
  }
}

setnames(master, c("XTY_A1", "XTY_A2"), c("A1", "A2"))
cols_base <- c("CHR", "SNP", "A1", "A2")
cols_all <- c(cols_base, sort(setdiff(names(master), cols_base)))
master <- master[, ..cols_all]
fwrite(master, file = final_out, sep = "\t", quote = FALSE, row.names = FALSE, na = "0")
EOF_RSCRIPT

# --- Run R Script ---
ORIGINAL_METHODS_STRING=$(IFS=,; echo "${ORIGINAL_METHODS_LIST[*]}")
Rscript "${R_SCRIPT_PATH}" \
    "${TARGET_TRAIT}" \
    "${TARGET_POPULATION_SHELL}" \
    "${CASE_DIR_NAME_CONCAT}" \
    "${CASE_TAG_CONCAT}" \
    "${REF_XTY_FILE}" \
    "${CONCAT_RESULTS_BASE_DIR}" \
    "${FINAL_OUTPUT_FILE_PATH}" \
    "${SOURCE_POP_LIST_FILE}" \
    "${ORIGINAL_METHODS_STRING}"

# --- Clean up ---
[ $? -eq 0 ] && rm -rf "${TEMP_DIR}"
