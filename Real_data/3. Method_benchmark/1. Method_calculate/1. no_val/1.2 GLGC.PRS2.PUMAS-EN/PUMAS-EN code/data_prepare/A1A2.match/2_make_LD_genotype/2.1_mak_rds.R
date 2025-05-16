library(bigsnpr)
library(fs)
library(purrr)

# --- Configuration ---
root_dir <- "/home/yd357/pi_paths/pi_zhao/GWAS_Subsample_JointPRS/v2_PUMAS-EN/PUMAS-EN/data/LD_ref/A1A2.match/2_make_LD_genotype"

# --- Locate All QCed Directories ---
qced_dirs <- dir_ls(root_dir, recurse = TRUE, type = "directory", regexp = "qced$")

# --- Conversion Function ---
process_qced <- function(qced_path) {
  bed_files <- dir_ls(qced_path, regexp = "\\.bed$")
  if (length(bed_files) == 0) return(NULL)
  
  bedfile <- bed_files[1]
  base_name <- path_ext_remove(path_file(bedfile))
  backingfile <- file.path(qced_path, base_name)
  
  if (file_exists(paste0(backingfile, ".rds"))) return(NULL)
  
  snp_readBed(bedfile, backingfile = backingfile)
}

# --- Convert All QCed PLINK Sets ---
walk(qced_dirs, process_qced)
