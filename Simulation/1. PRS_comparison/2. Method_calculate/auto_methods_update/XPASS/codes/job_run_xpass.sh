#!/bin/bash
#SBATCH --job-name=XPASS_array%j
#SBATCH --output=/home/yd357/pi_paths/pi_zhao/MIXPRS_Revise/Project1_15K_80K/code/4_methods/XPASS/logs/XPASS_array%j.log
#SBATCH --error=/home/yd357/pi_paths/pi_zhao/MIXPRS_Revise/Project1_15K_80K/code/4_methods/XPASS/logs/XPASS_array%j.log
#SBATCH --requeue
#SBATCH --mem=32G
#SBATCH -p scavenge,day,pi_zhao
#SBATCH -t 24:00:00
#SBATCH -c 1

# Load modules required by XPASS and its dependencies (ieugwasr::ld_clump)
module load R/4.3.0-foss-2020b
module load PLINK/1.9b_6.21-x86_64

# --- Paths ---
CODE_DIR="/home/yd357/pi_paths/pi_zhao/MIXPRS_Revise/Project1_15K_80K/code/4_methods/XPASS/codes"
INPUT_DIR="/home/yd357/pi_paths/pi_zhao/MIXPRS_Revise/Project1_15K_80K/data/sim_data/summary_data/discover_validate/XPASS"
RESULT_DIR="/home/yd357/pi_paths/pi_zhao/MIXPRS_Revise/Project1_15K_80K/result/add_benchmark/25.10.4/sim_data/XPASS"
# 1000G reference genotype paths from example
REF_DIR_BASE="/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/1000g_phase3_data/geno_data"

# --- Read Parameters ---
IFS=',' read -r SIM_I P SAMPLE_SIZE_2 POP2 < <(sed -n "${SLURM_ARRAY_TASK_ID}p" ${CODE_DIR}/xpass_params.txt)
RHOG="0.8"
POP1="EUR"
SAMPLE_SIZE_1="80K" # EUR label

# --- Define I/O files ---
# XPASS: file_z1 is target (pop2), file_z2 is auxiliary (pop1)
SS2_FILE_TARGET="${INPUT_DIR}/${POP2}_sim${SIM_I}_p${P}_rho${RHOG}_${SAMPLE_SIZE_2}_subcol.txt" # Target (POP2)
SS1_FILE_AUX="${INPUT_DIR}/${POP1}_sim${SIM_I}_p${P}_rho${RHOG}_${SAMPLE_SIZE_1}_subcol.txt" # Aux (EUR)

# XPASS: file_ref1 is target (pop2), file_ref2 is auxiliary (pop1)
REF_POP2_BFILE_TARGET="${REF_DIR_BASE}/${POP2}" # Target (POP2)
REF_POP1_BFILE_AUX="${REF_DIR_BASE}/${POP1}" # Aux (EUR)

OUT_FILE="${RESULT_DIR}/sim${SIM_I}_p${P}_rho${RHOG}_${SAMPLE_SIZE_2}_XPASS_EUR_${POP2}_beta_${POP2}.txt"

mkdir -p "${RESULT_DIR}"

echo "========================================================"
echo "RUNNING XPASS with parameters:"
echo "  Simulation ID : ${SIM_I}"
echo "  P-value       : ${P}"
echo "  Sample Scen.  : ${SAMPLE_SIZE_2} (for ${POP2})"
echo "  Target Pop    : ${POP2} (file: ${SS2_FILE_TARGET})"
echo "  Auxiliary Pop : ${POP1} (file: ${SS1_FILE_AUX})"
echo "  Target Ref    : ${REF_POP2_BFILE_TARGET}"
echo "  Auxiliary Ref : ${REF_POP1_BFILE_AUX}"
echo "  Output File   : ${OUT_FILE}"
echo "========================================================"

# Run R script via HEREDOC
# We pass paths and the pop2 name as command-line arguments to R
Rscript --vanilla - "$SS2_FILE_TARGET" "$SS1_FILE_AUX" "$REF_POP2_BFILE_TARGET" "$REF_POP1_BFILE_AUX" "$OUT_FILE" "$POP2" << 'EOF'
library(XPASS)
library(ieugwasr)
library(data.table)

# Get command line arguments
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 6) {
  stop("Usage: Rscript script.R <file_z1_target> <file_z2_aux> <file_ref1_target> <file_ref2_aux> <out_file> <pop2_name>", call.=FALSE)
}

file_z1_target_path <- args[1]
file_z2_aux_path <- args[2]
file_ref1_target_path <- args[3]
file_ref2_aux_path <- args[4]
out_file_path <- args[5]
pop2_name <- args[6]

# --- Clumping for Target Population (Pop1 in XPASS function, Pop2 in our var names) ---
message("Reading target summary stats: ", file_z1_target_path)
z_target <- fread(file_z1_target_path)
# Assume file has SNP and Z columns, per example
pval_target <- data.frame(rsid=z_target$SNP, pval=2*pnorm(abs(z_target$Z), lower.tail=FALSE))

snps_fe1 <- NULL
if(length(which(pval_target$pval < 1e-10)) > 0){
  message("Clumping target population (", pop2_name, ") SNPs...")
  tryCatch({
    clp_target <- ld_clump(pval_target, clump_kb=1000, clump_r2=0.1, clump_p=1e-10,
                           bfile=file_ref1_target_path, plink_bin="plink")
    snps_fe1 <- clp_target$rsid
    message("Found ", length(snps_fe1), " SNPs to use for fixed effects in target pop.")
  }, error = function(e) {
    message("LD clumping failed for target pop. Proceeding with snps_fe1=NULL. Error: ", e$message)
  })
} else {
  message("No SNPs below p=1e-10 in target pop, snps_fe1=NULL")
}

# --- Clumping for Auxiliary Population (Pop2 in XPASS function, Pop1 in our var names) ---
message("Reading auxiliary summary stats: ", file_z2_aux_path)
z_aux <- fread(file_z2_aux_path)
pval_aux <- data.frame(rsid=z_aux$SNP, pval=2*pnorm(abs(z_aux$Z), lower.tail=FALSE))

snps_fe2 <- NULL
if(length(which(pval_aux$pval < 1e-10)) > 0){
  message("Clumping auxiliary population (EUR) SNPs...")
   tryCatch({
    clp_aux <- ld_clump(pval_aux, clump_kb=1000, clump_r2=0.1, clump_p=1e-10,
                          bfile=file_ref2_aux_path, plink_bin="plink")
    snps_fe2 <- clp_aux$rsid
    message("Found ", length(snps_fe2), " SNPs to use for fixed effects in auxiliary pop.")
  }, error = function(e) {
    message("LD clumping failed for auxiliary pop. Proceeding with snps_fe2=NULL. Error: ", e$message)
  })
} else {
  message("No SNPs below p=1e-10 in auxiliary pop, snps_fe2=NULL")
}

# --- Run XPASS ---
message("Running XPASS...")
# file_z1 = target (pop2)
# file_z2 = auxiliary (pop1)
# file_ref1 = target (pop2)
# file_ref2 = auxiliary (pop1)
# snps_fe1 = target (pop2)
# snps_fe2 = auxiliary (pop1)
post_beta <- XPASS(file_z1 = file_z1_target_path, file_z2 = file_z2_aux_path,
                   file_ref1 = file_ref1_target_path, file_ref2 = file_ref2_aux_path,
                   snps_fe1 = snps_fe1, snps_fe2 = snps_fe2)

# --- Format and Write Output ---
message("Formatting and writing output to: ", out_file_path)
# mu_XPASS1 is the posterior mean for the target population (pop1 in XPASS terms)
final_beta <- post_beta$mu[, c("SNP", "A1", "mu_XPASS1")]
colnames(final_beta) <- c("rsID", "A1", pop2_name)

write.table(final_beta, file=out_file_path, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

message("XPASS complete.")
EOF

echo "R script finished for ${SIM_I}, ${P}, ${SAMPLE_SIZE_2}, ${POP2}"