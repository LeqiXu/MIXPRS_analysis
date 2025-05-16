#!/usr/bin/env Rscript

suppressMessages({
  library(optparse)
  library(bigsnpr)
  library(bigreadr)
  library(data.table)
  library(tidyverse)
})

# --- Parse Arguments ---
option_list <- list(
  make_option(c("-g", "--gwas"), type = "character"),
  make_option(c("-l", "--ld_ref"), type = "character"),
  make_option(c("-o", "--out"), type = "character"),
  make_option(c("--chr"), type = "integer"),
  make_option(c("--h2"), type = "numeric"),
  make_option(c("--total_snps"), type = "integer"),
  make_option(c("--ancestry"), type = "character", default = "unknown"),
  make_option(c("--cores"), type = "integer", default = 1),
  make_option(c("--mode"), type = "character"),
  make_option(c("--p"), type = "numeric", default = NA),
  make_option(c("--h2_ratio"), type = "numeric", default = NA),
  make_option(c("--sparse"), type = "character", default = NA),
  make_option(c("--shrink_corr"), type = "numeric", default = 0.95)
)
opt <- parse_args(OptionParser(option_list = option_list))
NCORES <- opt$cores
setDTthreads(NCORES)
Sys.setenv(OPENBLAS_NUM_THREADS = NCORES, MKL_NUM_THREADS = NCORES, OMP_NUM_THREADS = NCORES)

# --- Load LD Reference ---
rds_file <- paste0(opt$ld_ref, ".rds")
bed_file <- paste0(opt$ld_ref, ".bed")
if (!file.exists(rds_file)) {
  snp_readBed(bed_file, backingfile = sub("\\.rds$", "", rds_file))
}
ld_ref_obj <- snp_attach(rds_file)
G <- ld_ref_obj$genotypes
map_ld_ref <- ld_ref_obj$map[, c("chromosome", "marker.ID", "physical.pos", "allele1", "allele2")]
colnames(map_ld_ref) <- c("chr", "rsid", "pos", "a1", "a0")
map_ld_ref_chr <- map_ld_ref[map_ld_ref$chr == opt$chr, ]

# --- Load and Format GWAS ---
summstats <- fread2(opt$gwas, col.names = c("chr", "rsid", "pos", "a1", "a0", "maf", "beta", "beta_se", "n_eff")) %>%
  mutate(across(c(chr, pos), as.integer),
         across(c(maf, beta, beta_se, n_eff), as.numeric)) %>%
  filter(!is.na(chr) & !is.na(pos) & !is.na(beta) & !is.na(beta_se) & !is.na(n_eff) & beta_se > 0)

# --- Match SNPs ---
info_snp <- snp_match(summstats, map_ld_ref_chr, join_by_pos = TRUE, match.min.prop = 0.05)
if (nrow(info_snp) == 0) {
  file.create(opt$out); quit(save = "no", status = 0)
}

# --- Heritability Adjustment ---
snp_ratio <- nrow(info_snp) / opt$total_snps
h2_est <- opt$h2 * snp_ratio

# --- SD Filtering ---
ind_col_ref <- info_snp$`_NUM_ID_`
maf_ld_ref <- snp_MAF(G, ind.col = ind_col_ref, ncores = NCORES)
sd_ld_ref <- sqrt(2 * maf_ld_ref * (1 - maf_ld_ref))
sd_gwas <- sqrt(info_snp$maf * 2 * (1 - info_snp$maf))
is_bad_sd <- abs(sd_gwas - sd_ld_ref) >= 0.05 | sd_gwas < 0.01 | sd_ld_ref < 0.01
df_beta <- info_snp[!is_bad_sd, c("beta", "beta_se", "n_eff", "_NUM_ID_")]
info_snp_filtered <- info_snp[!is_bad_sd, ]
if (nrow(df_beta) == 0) {
  file.create(opt$out); quit(save = "no", status = 0)
}

# --- Compute LD Matrix ---
ind_chr <- df_beta$`_NUM_ID_`
corr_chr <- snp_cor(G, ind.col = ind_chr, size = 6000, infos.pos = info_snp_filtered$pos, ncores = NCORES)
if (all(is.na(corr_chr))) {
  file.create(opt$out); quit(save = "no", status = 0)
}
corr_sp <- as_SFBM(corr_chr, compact = TRUE)

# --- Run LDPred2 ---
final_results_df <- NULL

if (h2_est < 1e-4) {
  final_results_df <- data.frame(
    rsid = info_snp_filtered$rsid,
    a1 = info_snp_filtered$a1,
    beta = rep(0, nrow(info_snp_filtered))
  )
} else if (opt$mode == "grid") {
  h2_target <- max(1e-5, min(h2_est, opt$h2_ratio * h2_est))
  grid_params <- data.frame(p = opt$p, h2 = h2_target, sparse = opt$sparse == "T")
  beta_grid <- snp_ldpred2_grid(corr_sp, df_beta, grid_params, ncores = NCORES)[, 1]
  beta_grid[is.na(beta_grid) | abs(beta_grid) >= 1] <- 0
  final_results_df <- data.frame(
    rsid = info_snp_filtered$rsid,
    a1 = info_snp_filtered$a1,
    beta = beta_grid
  )
} else if (opt$mode == "auto") {
  multi_auto <- snp_ldpred2_auto(
    corr = corr_sp,
    df_beta = df_beta,
    h2_init = h2_est,
    vec_p_init = seq_log(1e-4, 0.9, length.out = NCORES),
    allow_jump_sign = FALSE,
    shrink_corr = opt$shrink_corr,
    ncores = NCORES
  )
  beta_auto <- multi_auto[[1]]$beta_est
  beta_auto[is.na(beta_auto)] <- 0
  final_results_df <- data.frame(
    rsid = info_snp_filtered$rsid,
    a1 = info_snp_filtered$a1,
    beta = beta_auto
  )
}

# --- Write Output ---
if (!is.null(final_results_df)) {
  dir.create(dirname(opt$out), showWarnings = FALSE, recursive = TRUE)
  fwrite(final_results_df, opt$out, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE, na = "NA")
} else {
  file.create(opt$out)
}
