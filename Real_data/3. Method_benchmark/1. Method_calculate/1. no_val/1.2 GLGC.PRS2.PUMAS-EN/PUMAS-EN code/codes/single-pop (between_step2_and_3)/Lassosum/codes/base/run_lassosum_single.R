#!/usr/bin/env Rscript
suppressMessages({
  library(optparse)
  library(data.table)
  library(lassosum)
  library(methods)
})

# --- Options ---
option_list <- list(
  make_option(c("-g", "--gwas_file"), type="character"),
  make_option(c("-n", "--sample_size"), type="integer"),
  make_option(c("-r", "--ref_bfile"), type="character"),
  make_option(c("-l", "--ld_blocks"), type="character"),
  make_option(c("-s", "--s_param"), type="character", default="0.2,0.5,0.9"),
  make_option(c("--lambda"), type="character", default="0.005,0.01"),
  make_option(c("-o", "--out_prefix"), type="character")
)
opt <- parse_args(OptionParser(option_list=option_list))

# --- Parameters ---
s_params <- as.numeric(strsplit(opt$s_param, ",")[[1]])
lambda_params <- as.numeric(strsplit(opt$lambda, ",")[[1]])

# --- Load GWAS ---
summstats <- fread(opt$gwas_file)
min_p <- min(summstats$P[summstats$P > 0], na.rm = TRUE)
summstats$P[summstats$P == 0] <- min_p
cor <- p2cor(p = summstats$P, n = opt$sample_size, sign = summstats$BETA)

# --- Load LD Blocks ---
LDblocks <- fread(opt$ld_blocks)
ld_chr <- unique(summstats$CHR)[1]
LDblocks <- LDblocks[LDblocks$chr == paste0("chr", ld_chr), ]

# --- Run Lassosum ---
out <- lassosum.pipeline(
  cor = cor,
  chr = summstats$CHR,
  pos = NULL,
  snp = summstats$SNP,
  A1 = summstats$A1,
  A2 = summstats$A2,
  ref.bfile = opt$ref_bfile,
  test.bfile = opt$ref_bfile,
  LDblocks = LDblocks,
  lambda = lambda_params,
  s = s_params,
  destandardize = TRUE,
  cluster = NULL
)

# --- Save Output ---
beta_matrix <- as.data.frame(out$beta)
snp_info <- out$sumstats[, c("snp", "A1")]
output_data <- data.frame(
  SNP = snp_info$snp,
  A1 = snp_info$A1,
  BETA = beta_matrix[, 1]
)

s_tag <- paste0("sX", opt$s_param)
lambda_tag <- if (opt$lambda == "0.005") {
  "lambda_0p005"
} else if (opt$lambda == "0.01") {
  "lambda_0p01"
} else {
  paste0("lambda_unexpected_", gsub("\\.", "p", opt$lambda))
}

output_file <- paste0(opt$out_prefix, "_", s_tag, "_", lambda_tag, ".lassosum.betas.txt")
fwrite(output_data, output_file, sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
