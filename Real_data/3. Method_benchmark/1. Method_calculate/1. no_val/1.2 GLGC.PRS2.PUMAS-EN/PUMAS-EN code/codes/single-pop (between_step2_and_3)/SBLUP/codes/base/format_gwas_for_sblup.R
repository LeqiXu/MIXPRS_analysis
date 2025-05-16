#!/usr/bin/env Rscript
suppressMessages(library(data.table))
options(stringsAsFactors=FALSE)

# --- Argument Parsing ---
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 4) {
  stop("Usage: Rscript format_gwas_for_sblup.R <input_gwas> <output_ma> <chromosome> <output_snp_count_m>", call. = FALSE)
}
gwas_raw_file <- args[1]
gwas_ma_file <- args[2]
target_chr <- as.integer(args[3])
snp_count_file <- args[4]

# --- File Reading ---
gwas_raw <- fread(gwas_raw_file, header = TRUE, sep = "auto", na.strings = c("NA", "NaN", "", "."))

# --- Column Verification ---
expected_cols <- c("CHR", "SNP", "A1", "A2", "MAF", "BETA", "SE", "P", "N")
if (!all(expected_cols %in% colnames(gwas_raw))) {
  stop("Missing required columns: ", paste(setdiff(expected_cols, colnames(gwas_raw)), collapse=", "), call. = FALSE)
}

# --- Filtering by Chromosome ---
gcta_df <- gwas_raw[gwas_raw$CHR == target_chr, ]

# --- Data Frame Creation ---
if (nrow(gcta_df) > 0) {
  gcta_df <- data.frame(
    SNP = gcta_df$SNP,
    A1 = gcta_df$A1,
    A2 = gcta_df$A2,
    freq = as.numeric(gcta_df$MAF),
    b = as.numeric(gcta_df$BETA),
    se = as.numeric(gcta_df$SE),
    p = as.numeric(gcta_df$P),
    N = as.numeric(gcta_df$N),
    stringsAsFactors = FALSE
  )
  gcta_df <- gcta_df[complete.cases(gcta_df[, c("SNP", "A1", "A2", "freq", "b", "se", "p", "N")]), ]
} else {
  gcta_df <- data.frame(SNP=character(), A1=character(), A2=character(), freq=numeric(),
                        b=numeric(), se=numeric(), p=numeric(), N=numeric())
}

# --- File Writing ---
fwrite(gcta_df, gwas_ma_file, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE, na = "NA")
writeLines(as.character(nrow(gcta_df)), con = snp_count_file)
