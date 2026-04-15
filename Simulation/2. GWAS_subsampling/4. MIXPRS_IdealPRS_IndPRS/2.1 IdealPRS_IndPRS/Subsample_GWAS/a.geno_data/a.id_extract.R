## 0. Intersect with hapmaps
## hm3 SNP: 1203063
library(data.table)

snp_info = fread("/gpfs/gibbs/pi/zhao/xz674/data/sim_data/geno_data/info/snp_infor_mega+hm3")
snp_hm3 = fread("/gpfs/gibbs/pi/zhao/gz222/1000g_phase3/Hapmap3_snp/hm3_noMHC.snplist", header = F)
snp_info_hm3 = intersect(snp_info$rs_id,snp_hm3$V1)
SNP_info_hm3_SNP = snp_info$SNP[which(snp_info$rs_id %in% snp_info_hm3)]
SNP_info_hm3_rs_id = snp_info$rs_id[which(snp_info$rs_id %in% snp_info_hm3)]

SNP_map = data.table(SNP = SNP_info_hm3_SNP, rs_id = SNP_info_hm3_rs_id)
print(paste0("hm3 snp size: ", length(SNP_info_hm3_SNP)))

write.table(SNP_info_hm3_SNP, "/gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/info/snp_infor_hm3.snplist", row.names=F, col.names=F, quote=F, append=F)
write.table(SNP_info_hm3_rs_id, "/gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/info/snp_rs_id_infor_hm3.snplist", row.names=F, col.names=F, quote=F, append=F)
write.table(SNP_map, "/gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/info/snp_map_hm3.snplist", row.names=F, col.names=F, quote=F, append=F, sep = "\t")


## 1. Extract discover and validate people from simulated data
n_fold = 4

# test
test_id = data.table(FID = c(110001:120000), IID = c(110001:120000))
write.table(test_id, "/gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/info/test_id.tsv", row.names=F, col.names=T, quote = F, sep = "\t")

# whole
discover_validate_20K_id = data.table(FID = c(1:20000), IID = c(1:20000))
discover_validate_100K_id = data.table(FID = c(1:100000), IID = c(1:100000))

write.table(discover_validate_20K_id, "/gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/info/discover_validate_20K_id_all.tsv", row.names=F, col.names=T, quote = F, sep = "\t")
write.table(discover_validate_100K_id, "/gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/info/discover_validate_100K_id_all.tsv", row.names=F, col.names=T, quote = F, sep = "\t")

# cross-valid split
fold_20K_id <- split(1:nrow(discover_validate_20K_id), cut(1:nrow(discover_validate_20K_id), breaks = n_fold, labels = FALSE))
fold_100K_id <- split(1:nrow(discover_validate_100K_id), cut(1:nrow(discover_validate_100K_id), breaks = n_fold, labels = FALSE))

for (ff in 1:n_fold){
  validate_20K_id = discover_validate_20K_id[fold_20K_id[[ff]], ]
  validate_100K_id = discover_validate_100K_id[fold_100K_id[[ff]], ]
  discover_20K_id = discover_validate_20K_id[-fold_20K_id[[ff]], ]
  discover_100K_id = discover_validate_100K_id[-fold_100K_id[[ff]], ]
  
  write.table(validate_20K_id, paste0("/gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/info/validate_20K_id_fold",ff,".tsv"), row.names=F, col.names=T, quote = F, sep = "\t")
  write.table(validate_100K_id, paste0("/gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/info/validate_100K_id_fold",ff,".tsv"), row.names=F, col.names=T, quote = F, sep = "\t")
  write.table(discover_20K_id, paste0("/gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/info/discover_20K_id_fold",ff,".tsv"), row.names=F, col.names=T, quote = F, sep = "\t")
  write.table(discover_100K_id, paste0("/gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/info/discover_100K_id_fold",ff,".tsv"), row.names=F, col.names=T, quote = F, sep = "\t")
}


## 2. Extract discover and validate people from 1KG data

# For each population, obtain their id
for (pop in c("EUR","EAS","AFR","SAS","AMR")){
  pop_1kg_id = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/1000g_phase3_data/ancestry_info/",pop,"_id.tsv"))
  
  # Randomly sample the ids for 1kg1, and the others for 1kg2
  idx_1kg1 = sort(sample(1:nrow(pop_1kg_id), ceiling(nrow(pop_1kg_id) / 2), replace = FALSE))
  pop_1kg1_id = pop_1kg_id[idx_1kg1, ]
  pop_1kg2_id = pop_1kg_id[- idx_1kg1, ]
  
  # whole
  write.table(pop_1kg_id, paste0("/gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/info/",pop,"_discover_validate_1kg_id_all.tsv"), row.names=F, col.names=T, quote = F, sep = "\t")
  write.table(pop_1kg1_id, paste0("/gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/info/",pop,"_discover_validate_1kg1_id_all.tsv"), row.names=F, col.names=T, quote = F, sep = "\t")
  write.table(pop_1kg2_id, paste0("/gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/info/",pop,"_discover_validate_1kg2_id_all.tsv"), row.names=F, col.names=T, quote = F, sep = "\t")
  
  # cross-valid split
  fold_1kg_id <- split(1:nrow(pop_1kg_id), cut(1:nrow(pop_1kg_id), breaks = n_fold, labels = FALSE))
  fold_1kg1_id <- split(1:nrow(pop_1kg1_id), cut(1:nrow(pop_1kg1_id), breaks = n_fold, labels = FALSE))
  fold_1kg2_id <- split(1:nrow(pop_1kg2_id), cut(1:nrow(pop_1kg2_id), breaks = n_fold, labels = FALSE))
  
  for (ff in 1:n_fold){
    validate_1kg_id = pop_1kg_id[fold_1kg_id[[ff]], ]
    validate_1kg1_id = pop_1kg1_id[fold_1kg1_id[[ff]], ]
    validate_1kg2_id = pop_1kg2_id[fold_1kg2_id[[ff]], ]
    discover_1kg_id = pop_1kg_id[-fold_1kg_id[[ff]], ]
    discover_1kg1_id = pop_1kg1_id[-fold_1kg1_id[[ff]], ]
    discover_1kg2_id = pop_1kg2_id[-fold_1kg2_id[[ff]], ]
    
    write.table(validate_1kg_id, paste0("/gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/info/",pop,"_validate_1kg_id_fold",ff,".tsv"), row.names=F, col.names=T, quote = F, sep = "\t")
    write.table(validate_1kg1_id, paste0("/gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/info/",pop,"_validate_1kg1_id_fold",ff,".tsv"), row.names=F, col.names=T, quote = F, sep = "\t")
    write.table(validate_1kg2_id, paste0("/gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/info/",pop,"_validate_1kg2_id_fold",ff,".tsv"), row.names=F, col.names=T, quote = F, sep = "\t")
    write.table(discover_1kg_id, paste0("/gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/info/",pop,"_discover_1kg_id_fold",ff,".tsv"), row.names=F, col.names=T, quote = F, sep = "\t")
    write.table(discover_1kg1_id, paste0("/gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/info/",pop,"_discover_1kg1_id_fold",ff,".tsv"), row.names=F, col.names=T, quote = F, sep = "\t")
    write.table(discover_1kg2_id, paste0("/gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/info/",pop,"_discover_1kg2_id_fold",ff,".tsv"), row.names=F, col.names=T, quote = F, sep = "\t")
  }
}


## 3. Extract discover and validate people from UK Biobank data

# For each population, obtain their id
for (pop in c("EUR")){
  pop_ukbb_id = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/ancestry_info/",pop,"_id.tsv"))

  # whole
  write.table(pop_ukbb_id, paste0("/gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/info/",pop,"_discover_validate_ukbb_id_all.tsv"), row.names=F, col.names=T, quote = F, sep = "\t")

  # cross-valid split
  fold_ukbb_id <- split(1:nrow(pop_ukbb_id), cut(1:nrow(pop_ukbb_id), breaks = n_fold, labels = FALSE))
  
  for (ff in 1:n_fold){
    validate_ukbb_id = pop_ukbb_id[fold_ukbb_id[[ff]], ]
    discover_ukbb_id = pop_ukbb_id[-fold_ukbb_id[[ff]], ]
    
    write.table(validate_ukbb_id, paste0("/gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/info/",pop,"_validate_ukbb_id_fold",ff,".tsv"), row.names=F, col.names=T, quote = F, sep = "\t")
    write.table(discover_ukbb_id, paste0("/gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/info/",pop,"_discover_ukbb_id_fold",ff,".tsv"), row.names=F, col.names=T, quote = F, sep = "\t")
  }
}

