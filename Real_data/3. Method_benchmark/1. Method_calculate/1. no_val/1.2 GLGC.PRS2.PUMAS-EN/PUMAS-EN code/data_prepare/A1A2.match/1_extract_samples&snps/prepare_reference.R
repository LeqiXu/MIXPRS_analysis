library(data.table)
options(stringsAsFactors = FALSE)

# --- Population Config ---
pop_codes <- c("EUR", "EAS", "AFR", "AMR", "SAS")
pop_names <- tolower(pop_codes)

# --- Study Config ---
studies <- list(
  GLGC = list(pops = 1:5),
  PAGE = list(pops = 1:3),
  BBJ  = list(pops = 1:2)
)

base_output_dir <- "/home/yd357/pi_paths/pi_zhao/GWAS_Subsample_JointPRS/v2_PUMAS-EN/PUMAS-EN/data/LD_ref/A1A2.match/1_extract_samples&snps"
dir.create(base_output_dir, recursive = TRUE, showWarnings = FALSE)

# --- Population Processing Function ---
process_population <- function(study_name, pop_index) {
  pop <- pop_codes[pop_index]
  pop_name <- pop_names[pop_index]
  
  study_dir <- file.path(base_output_dir, study_name)
  tmp_dir <- file.path(study_dir, "tmp")
  
  geno_path <- file.path("/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/1000g_phase3_data/geno_data", pop)
  ind <- fread(paste0(geno_path, ".fam"), header = FALSE)
  
  sample_size <- min(nrow(ind), 500)
  tmp.ind <- ind[sample(nrow(ind), size = sample_size), .(V1, V2)]
  fwrite(tmp.ind, file.path(tmp_dir, paste0(pop, "_LD500.samples.txt")),
         col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
  
  cs_file <- paste0("/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ref_data/PRScsx/1kg/ldblk_1kg_", pop_name, "/snpinfo_1kg_hm3")
  CS <- fread(cs_file)
  LD <- fread(paste0(geno_path, ".bim"))
  
  snp <- intersect(LD$V2, CS$SNP)
  LD <- LD[match(snp, LD$V2)]
  CS <- CS[match(snp, CS$SNP)]
  keep <- (LD$V5 == CS$A1 & LD$V6 == CS$A2) | (LD$V6 == CS$A1 & LD$V5 == CS$A2)
  CS <- CS[keep]
  
  fwrite(CS[, .(SNP, A1)], file.path(study_dir, paste0("PRScs_", pop_name, "_1kg_A1.txt")),
         col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
  fwrite(CS[, .(SNP)], file.path(study_dir, paste0("PRScs_", pop_name, "_1kg.snplist")),
         col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
}

# --- Process Each Study ---
for (study_name in names(studies)) {
  study_dir <- file.path(base_output_dir, study_name)
  tmp_dir <- file.path(study_dir, "tmp")
  dir.create(study_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
  
  for (pop_index in studies[[study_name]]$pops) {
    process_population(study_name, pop_index)
  }
}
