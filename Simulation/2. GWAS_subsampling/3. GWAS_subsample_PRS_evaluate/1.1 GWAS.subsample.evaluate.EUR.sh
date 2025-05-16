############################## Step0: GWAS analysis function ##############################
# Check the complement
complement <- function(a_vec) {
  # Ensure uppercase
  a_vec <- toupper(a_vec)
  # Create a lookup
  complement_map <- c(A = "T", T = "A", C = "G", G = "C")
  out <- complement_map[a_vec]
  # For any that are NA in the map, keep the original or handle as needed
  out[is.na(out)] <- a_vec[is.na(out)]
  return(out)
}

# Flip function
flip_merge_effect <- function(gwas_df, true_effect_df) {
  # Merge with true effect
  merged_df <- merge(gwas_df[, c("SNP","A1","BETA_STD","N")], true_effect_df[, c("true_effect_SNP","true_effect_A1","true_effect_BETA_STD")], by.x = "SNP", by.y = "true_effect_SNP", all.x = TRUE)
  
  # Check direct match OR complement match
  same_allele <- (merged_df$A1 == merged_df$true_effect_A1)
  comp_allele <- (merged_df$A1 == complement(merged_df$true_effect_A1))
  match_allele <- same_allele | comp_allele
  
  # Flip sign if no match
  merged_df$true_effect_BETA_STD <- ifelse(match_allele, merged_df$true_effect_BETA_STD, -merged_df$true_effect_BETA_STD)
  
  # Set effect = 0 for SNPs not in 'true_effect'
  merged_df$true_effect_BETA_STD[ !(merged_df$SNP %in% true_effect_df$true_effect_SNP) ] <- 0
  
  return(merged_df)
}

# LD block eigen decomp function
make_svd_averaged <- function(mat) {
  # Perform SVD
  svd_decomp <- svd(mat)
  V <- svd_decomp$v
  s <- svd_decomp$d

  # Reconstruct symmetric approximation
  h <- V %*% diag(s) %*% t(V)

  # Average with original
  mat_sym <- (mat + h) / 2
  return(mat_sym)
}

eigen_decomp_function = function(ldblk){

  eigen_decomp = eigen(ldblk, symmetric = TRUE)  # Use symmetric = TRUE for efficiency
  Q = eigen_decomp$vectors
  Lambda = eigen_decomp$values
    
  # Compute sqrt and inverse sqrt with thresholding
  threshold = 1e-10
  Lambda_sqrt = pmax(Lambda, threshold) ^ 0.5
  ldblk_half = Q %*% diag(Lambda_sqrt) %*% t(Q)

  return(ldblk_half)
}

eigen_decomp_inv_function = function(ldblk){

  eigen_decomp = eigen(ldblk, symmetric = TRUE)  # Use symmetric = TRUE for efficiency
  Q = eigen_decomp$vectors
  Lambda = eigen_decomp$values
    
  ### Compute sqrt and inverse sqrt with thresholding
  threshold = 1e-10
  Lambda_inv_sqrt = pmax(Lambda, threshold) ^ -0.5
  ldblk_inv_half = Q %*% diag(Lambda_inv_sqrt) %*% t(Q)

  return(ldblk_inv_half)
}

# residual calculate function
residual_function_identity = function(ldblk_sim_ref, ldblk_sim, snplist_sim, ldblk_1kg_ref, ldblk_1kg, snplist_1kg, snplist_inter, subsample_GWAS) {
  
  # Subset subsample_GWAS and ldblk
  subsample_GWAS_sub = subsample_GWAS[match(snplist_inter, subsample_GWAS$SNP), ]
  snplist_sim_ld_idx = match(snplist_inter, snplist_sim)  
  snplist_1kg_ld_idx = match(snplist_inter, snplist_1kg)  

  ldblk_sim = ldblk_sim[snplist_sim_ld_idx, snplist_sim_ld_idx]
  ldblk_1kg = ldblk_1kg[snplist_1kg_ld_idx, snplist_1kg_ld_idx]

  ldblk_sim_ref_sub = ldblk_sim_ref[match(snplist_inter, ldblk_sim_ref$ref_sim_SNP),]
  ldblk_1kg_ref_sub = ldblk_1kg_ref[match(snplist_inter, ldblk_1kg_ref$ref_1kg_SNP),]

  # ldblk_sim and ldblk_1kg
  ldblk_sim_same_or_complement <- (subsample_GWAS_sub$A1 == ldblk_sim_ref_sub$ref_sim_A1) | (subsample_GWAS_sub$A1 == complement(ldblk_sim_ref_sub$ref_sim_A1))
  ldblk_1kg_same_or_complement <- (subsample_GWAS_sub$A1 == ldblk_1kg_ref_sub$ref_1kg_A1) | (subsample_GWAS_sub$A1 == complement(ldblk_1kg_ref_sub$ref_1kg_A1))

  ldblk_sim_mismatch_id = which(!ldblk_sim_same_or_complement)
  ldblk_1kg_mismatch_id = which(!ldblk_1kg_same_or_complement)

  if(length(ldblk_sim_mismatch_id) > 0 & length(snplist_inter) > 1){
    ldblk_sim[ldblk_sim_mismatch_id,] = -ldblk_sim[ldblk_sim_mismatch_id,]
    ldblk_sim[,ldblk_sim_mismatch_id] = -ldblk_sim[,ldblk_sim_mismatch_id]
  }
  if(length(ldblk_1kg_mismatch_id) > 0 & length(snplist_inter) > 1){
    ldblk_1kg[ldblk_1kg_mismatch_id,] = -ldblk_1kg[ldblk_1kg_mismatch_id,]
    ldblk_1kg[,ldblk_1kg_mismatch_id] = -ldblk_1kg[,ldblk_1kg_mismatch_id]
  }

  current_N = median(subsample_GWAS_sub$N)

  ldblk_sim <- make_svd_averaged(ldblk_sim)
  ldblk_1kg <- make_svd_averaged(ldblk_1kg)
  identity_matrix <- diag(nrow(ldblk_1kg))
  ldblk_current = current_N / GWAS_subsample_N * ldblk_sim + (GWAS_subsample_N - current_N) / GWAS_subsample_N * identity_matrix

  if (length(snplist_inter) > 1){
    # Eigenvalue decomposition    
    ldblk_inv_current_half = eigen_decomp_inv_function(ldblk_current)
    ldblk_sim_half = eigen_decomp_function(ldblk_sim)

    # Compute residuals directly
    subsample_GWAS_sub$Residual = as.vector(ldblk_inv_current_half %*% (subsample_GWAS_sub$BETA_STD - ldblk_sim_half %*% ldblk_sim_half %*% subsample_GWAS_sub$true_effect_BETA_STD))

  } else{
    # Compute residuals directly
    subsample_GWAS_sub$Residual = as.vector(1/sqrt(ldblk_current) * (subsample_GWAS_sub$BETA_STD - ldblk_sim * subsample_GWAS_sub$true_effect_BETA_STD))
    
  }

  subsample_GWAS_sub$Residual = sqrt(current_N) * subsample_GWAS_sub$Residual
  
  return(subsample_GWAS_sub)

}

residual_function_refLD = function(ldblk_sim_ref, ldblk_sim, snplist_sim, ldblk_1kg_ref, ldblk_1kg, snplist_1kg, snplist_inter, subsample_GWAS) {
  
  # Subset subsample_GWAS and ldblk
  subsample_GWAS_sub = subsample_GWAS[match(snplist_inter, subsample_GWAS$SNP), ]
  snplist_sim_ld_idx = match(snplist_inter, snplist_sim)  
  snplist_1kg_ld_idx = match(snplist_inter, snplist_1kg)  

  ldblk_sim = ldblk_sim[snplist_sim_ld_idx, snplist_sim_ld_idx]
  ldblk_1kg = ldblk_1kg[snplist_1kg_ld_idx, snplist_1kg_ld_idx]

  ldblk_sim_ref_sub = ldblk_sim_ref[match(snplist_inter, ldblk_sim_ref$ref_sim_SNP),]
  ldblk_1kg_ref_sub = ldblk_1kg_ref[match(snplist_inter, ldblk_1kg_ref$ref_1kg_SNP),]

  # ldblk_sim and ldblk_1kg
  ldblk_sim_same_or_complement <- (subsample_GWAS_sub$A1 == ldblk_sim_ref_sub$ref_sim_A1) | (subsample_GWAS_sub$A1 == complement(ldblk_sim_ref_sub$ref_sim_A1))
  ldblk_1kg_same_or_complement <- (subsample_GWAS_sub$A1 == ldblk_1kg_ref_sub$ref_1kg_A1) | (subsample_GWAS_sub$A1 == complement(ldblk_1kg_ref_sub$ref_1kg_A1))

  ldblk_sim_mismatch_id = which(!ldblk_sim_same_or_complement)
  ldblk_1kg_mismatch_id = which(!ldblk_1kg_same_or_complement)

  if(length(ldblk_sim_mismatch_id) > 0 & length(snplist_inter) > 1){
    ldblk_sim[ldblk_sim_mismatch_id,] = -ldblk_sim[ldblk_sim_mismatch_id,]
    ldblk_sim[,ldblk_sim_mismatch_id] = -ldblk_sim[,ldblk_sim_mismatch_id]
  }
  if(length(ldblk_1kg_mismatch_id) > 0 & length(snplist_inter) > 1){
    ldblk_1kg[ldblk_1kg_mismatch_id,] = -ldblk_1kg[ldblk_1kg_mismatch_id,]
    ldblk_1kg[,ldblk_1kg_mismatch_id] = -ldblk_1kg[,ldblk_1kg_mismatch_id]
  }

  current_N = median(subsample_GWAS_sub$N)

  ldblk_sim <- make_svd_averaged(ldblk_sim)
  ldblk_1kg <- make_svd_averaged(ldblk_1kg)
  ldblk_current = current_N / GWAS_subsample_N * ldblk_sim + (GWAS_subsample_N - current_N) / GWAS_subsample_N * ldblk_1kg

  if (length(snplist_inter) > 1){
    # Eigenvalue decomposition    
    ldblk_inv_current_half = eigen_decomp_inv_function(ldblk_current)
    ldblk_sim_half = eigen_decomp_function(ldblk_sim)

    # Compute residuals directly
    subsample_GWAS_sub$Residual = as.vector(ldblk_inv_current_half %*% (subsample_GWAS_sub$BETA_STD - ldblk_sim_half %*% ldblk_sim_half %*% subsample_GWAS_sub$true_effect_BETA_STD))

  } else{
    # Compute residuals directly
    subsample_GWAS_sub$Residual = as.vector(1/sqrt(ldblk_current) * (subsample_GWAS_sub$BETA_STD - ldblk_sim * subsample_GWAS_sub$true_effect_BETA_STD))
    
  }

  subsample_GWAS_sub$Residual = sqrt(current_N) * subsample_GWAS_sub$Residual
  
  return(subsample_GWAS_sub)

}

# GWAS residual analize function
GWAS_residual_analysis_function <- function(pop,true_effect,train_GWAS_full_snplist,tune_GWAS_full_snplist,train_GWAS_prune_snplist,tune_GWAS_prune_snplist,train_GWAS_prune_snplist_ind_approx,tune_GWAS_prune_snplist_ind_approx){

### Step1: Obtain standardized beta for subsampled GWAS
train_GWAS_full_snplist$BETA_STD = sign(train_GWAS_full_snplist$BETA) * abs(qnorm(train_GWAS_full_snplist$P / 2)) / sqrt(train_GWAS_full_snplist$N)
tune_GWAS_full_snplist$BETA_STD = sign(tune_GWAS_full_snplist$BETA) * abs(qnorm(tune_GWAS_full_snplist$P / 2)) / sqrt(tune_GWAS_full_snplist$N)
train_GWAS_prune_snplist$BETA_STD = sign(train_GWAS_prune_snplist$BETA) * abs(qnorm(train_GWAS_prune_snplist$P / 2)) / sqrt(train_GWAS_prune_snplist$N)
tune_GWAS_prune_snplist$BETA_STD = sign(tune_GWAS_prune_snplist$BETA) * abs(qnorm(tune_GWAS_prune_snplist$P / 2)) / sqrt(tune_GWAS_prune_snplist$N)
train_GWAS_prune_snplist_ind_approx$BETA_STD = sign(train_GWAS_prune_snplist_ind_approx$BETA) * abs(qnorm(train_GWAS_prune_snplist_ind_approx$P / 2)) / sqrt(train_GWAS_prune_snplist_ind_approx$N)
tune_GWAS_prune_snplist_ind_approx$BETA_STD = sign(tune_GWAS_prune_snplist_ind_approx$BETA) * abs(qnorm(tune_GWAS_prune_snplist_ind_approx$P / 2)) / sqrt(tune_GWAS_prune_snplist_ind_approx$N)

### Step2: Merge true effect with standardized beta
train_GWAS_full_snplist_true_effect  <- flip_merge_effect(train_GWAS_full_snplist,  true_effect)
tune_GWAS_full_snplist_true_effect   <- flip_merge_effect(tune_GWAS_full_snplist,   true_effect)
train_GWAS_prune_snplist_true_effect <- flip_merge_effect(train_GWAS_prune_snplist, true_effect)
tune_GWAS_prune_snplist_true_effect  <- flip_merge_effect(tune_GWAS_prune_snplist,  true_effect)
train_GWAS_prune_snplist_ind_approx_true_effect <- flip_merge_effect(train_GWAS_prune_snplist_ind_approx, true_effect)
tune_GWAS_prune_snplist_ind_approx_true_effect  <- flip_merge_effect(tune_GWAS_prune_snplist_ind_approx,  true_effect)

### Step3: Calculate residuals
# Obtain A1 for LD reference panel 
ldblk_sim_ref = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ref_data/PRScsx/ukbb/ldblk_ukbb_", tolower(pop),"/snpinfo_ukbb_hm3"))
ldblk_sim_ref = ldblk_sim_ref[,c("SNP","A1")]
colnames(ldblk_sim_ref) = c("ref_sim_SNP","ref_sim_A1")

ldblk_1kg_ref = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ref_data/PRScsx/1kg/ldblk_1kg_", tolower(pop), "/snpinfo_1kg_hm3"))
ldblk_1kg_ref = ldblk_1kg_ref[,c("SNP","A1")]
colnames(ldblk_1kg_ref) = c("ref_1kg_SNP","ref_1kg_A1")

# Calculate residual
combined_residual = list()

for (data_type in c("train_GWAS_full_snplist_true_effect","tune_GWAS_full_snplist_true_effect","train_GWAS_prune_snplist_true_effect","tune_GWAS_prune_snplist_true_effect","train_GWAS_prune_snplist_ind_approx_true_effect","tune_GWAS_prune_snplist_ind_approx_true_effect")){

combined_residual[[data_type]] = data.table()
subsample_GWAS = get(data_type)

for (chr in 1:22) {
  file_path_sim = paste0("/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ref_data/PRScsx/ukbb/ldblk_ukbb_", tolower(pop),"/ldblk_ukbb_chr", chr, ".hdf5")
  contents_sim = h5ls(file_path_sim)

  file_path_1kg = paste0("/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ref_data/PRScsx/1kg/ldblk_1kg_", tolower(pop),"/ldblk_1kg_chr", chr, ".hdf5")
  contents_1kg = h5ls(file_path_1kg)

  for (block_name in unique(contents_sim$group)) {
    if (block_name != "/") {

      snplist_sim = h5read(file_path_sim, paste0(block_name, "/snplist"))
      snplist_1kg = h5read(file_path_1kg, paste0(block_name, "/snplist"))

      snplist_inter = intersect(snplist_sim, subsample_GWAS$SNP)
      snplist_inter = intersect(snplist_1kg, snplist_inter)

      if(length(snplist_inter) > 0){

      ldblk_sim = h5read(file_path_sim, paste0(block_name, "/ldblk"))
      ldblk_1kg = h5read(file_path_1kg, paste0(block_name, "/ldblk"))

      if (data_type %in% c("train_GWAS_prune_snplist_ind_approx_true_effect","tune_GWAS_prune_snplist_ind_approx_true_effect")){
        sub_residual_table = residual_function_identity(ldblk_sim_ref, ldblk_sim, snplist_sim, ldblk_1kg_ref, ldblk_1kg, snplist_1kg, snplist_inter, subsample_GWAS)
      } else {
        sub_residual_table = residual_function_refLD(ldblk_sim_ref, ldblk_sim, snplist_sim, ldblk_1kg_ref, ldblk_1kg, snplist_1kg, snplist_inter, subsample_GWAS)
      }
      combined_residual[[data_type]] = rbind(combined_residual[[data_type]], sub_residual_table)

      }      
    }
  }
}

}

return(combined_residual)

}

############################## Step1: Read in the GWAS ##############################
library(data.table)
library(rhdf5)
library(Matrix)

h2 = 0.4
rhog = 0.8
pop = "EUR"
sample_size = "ukbb"

args <- commandArgs(trailingOnly = TRUE)
sim_i <- as.integer(args[1])  #  1 to 5
repeat_i <- as.integer(args[2])  # 1 to 4
p <- as.numeric(args[3]) # 0.1 0.01 0.001 5e-04

# True effect
true_effect = fread(paste0("/gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/pheno_data/",pop,"/discover_validate/",pop,"_sim",sim_i,"_h2",h2,"_p",p,"_rhog",rhog,"_",sample_size,"_all.par"))
colnames(true_effect) = c("true_effect_SNP","true_effect_A1","true_effect_A1_Frq","true_effect_BETA_STD")

# Subsampled GWAS with full snplist
train_GWAS_full_snplist = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/sim_data/summary_data/subsample_data/clean/",pop,"_sim",sim_i,"_h2",h2,"_p",p,"_rhog",rhog,"_",sample_size,"_full_snplist_",pop,"_train_GWAS_approxFALSE_ratio3.00_repeat",repeat_i,".txt"))
tune_GWAS_full_snplist = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/sim_data/summary_data/subsample_data/clean/",pop,"_sim",sim_i,"_h2",h2,"_p",p,"_rhog",rhog,"_",sample_size,"_full_snplist_",pop,"_tune_GWAS_approxFALSE_ratio3.00_repeat",repeat_i,".txt"))

# Subsample GWAS with prune snplist
i=1

train_GWAS_prune_snplist = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/sim_data/summary_data/subsample_data/clean/",pop,"_sim",sim_i,"_h2",h2,"_p",p,"_rhog",rhog,"_",sample_size,"_prune_snplist_",i,"_",pop,"_train_GWAS_approxFALSE_ratio3.00_repeat",repeat_i,".txt"))
tune_GWAS_prune_snplist = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/sim_data/summary_data/subsample_data/clean/",pop,"_sim",sim_i,"_h2",h2,"_p",p,"_rhog",rhog,"_",sample_size,"_prune_snplist_",i,"_",pop,"_tune_GWAS_approxFALSE_ratio3.00_repeat",repeat_i,".txt"))

train_GWAS_prune_snplist_ind_approx = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/sim_data/summary_data/subsample_data/clean/",pop,"_sim",sim_i,"_h2",h2,"_p",p,"_rhog",rhog,"_",sample_size,"_prune_snplist_",i,"_",pop,"_train_GWAS_approxTRUE_ratio3.00_repeat",repeat_i,".txt"))
tune_GWAS_prune_snplist_ind_approx = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/sim_data/summary_data/subsample_data/clean/",pop,"_sim",sim_i,"_h2",h2,"_p",p,"_rhog",rhog,"_",sample_size,"_prune_snplist_",i,"_",pop,"_tune_GWAS_approxTRUE_ratio3.00_repeat",repeat_i,".txt"))

# Subsample sample size
if(sample_size == "ukbb"){
  GWAS_subsample_N = 311600
}

combined_residual = GWAS_residual_analysis_function(pop,true_effect,train_GWAS_full_snplist,tune_GWAS_full_snplist,train_GWAS_prune_snplist,tune_GWAS_prune_snplist,train_GWAS_prune_snplist_ind_approx,tune_GWAS_prune_snplist_ind_approx)

combined_residual_full_snplist_train = combined_residual[["train_GWAS_full_snplist_true_effect"]]
combined_residual_full_snplist_tune = combined_residual[["tune_GWAS_full_snplist_true_effect"]]
combined_residual_prune_snplist_train = combined_residual[["train_GWAS_prune_snplist_true_effect"]]
combined_residual_prune_snplist_tune = combined_residual[["tune_GWAS_prune_snplist_true_effect"]]
combined_residual_prune_snplist_ind_approx_train = combined_residual[["train_GWAS_prune_snplist_ind_approx_true_effect"]]
combined_residual_prune_snplist_ind_approx_tune = combined_residual[["tune_GWAS_prune_snplist_ind_approx_true_effect"]]

write.table(combined_residual_full_snplist_train,paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/residual/",pop,"_sim",sim_i,"_h2",h2,"_p",p,"_rhog",rhog,"_",sample_size,"_full_snplist_",pop,"_train_GWAS_residual_approxFALSE_ratio3.00_repeat",repeat_i,".txt"),quote=F,sep='\t',row.names=F,col.names=T)
write.table(combined_residual_full_snplist_tune,paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/residual/",pop,"_sim",sim_i,"_h2",h2,"_p",p,"_rhog",rhog,"_",sample_size,"_full_snplist_",pop,"_tune_GWAS_residual_approxFALSE_ratio3.00_repeat",repeat_i,".txt"),quote=F,sep='\t',row.names=F,col.names=T)

write.table(combined_residual_prune_snplist_train,paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/residual/",pop,"_sim",sim_i,"_h2",h2,"_p",p,"_rhog",rhog,"_",sample_size,"_prune_snplist_",i,"_",pop,"_train_GWAS_residual_approxFALSE_ratio3.00_repeat",repeat_i,".txt"),quote=F,sep='\t',row.names=F,col.names=T)
write.table(combined_residual_prune_snplist_tune,paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/residual/",pop,"_sim",sim_i,"_h2",h2,"_p",p,"_rhog",rhog,"_",sample_size,"_prune_snplist_",i,"_",pop,"_tune_GWAS_residual_approxFALSE_ratio3.00_repeat",repeat_i,".txt"),quote=F,sep='\t',row.names=F,col.names=T)

write.table(combined_residual_prune_snplist_ind_approx_train,paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/residual/",pop,"_sim",sim_i,"_h2",h2,"_p",p,"_rhog",rhog,"_",sample_size,"_prune_snplist_",i,"_",pop,"_train_GWAS_residual_approxTRUE_ratio3.00_repeat",repeat_i,".txt"),quote=F,sep='\t',row.names=F,col.names=T)
write.table(combined_residual_prune_snplist_ind_approx_tune,paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/residual/",pop,"_sim",sim_i,"_h2",h2,"_p",p,"_rhog",rhog,"_",sample_size,"_prune_snplist_",i,"_",pop,"_tune_GWAS_residual_approxTRUE_ratio3.00_repeat",repeat_i,".txt"),quote=F,sep='\t',row.names=F,col.names=T)


vim EUR_residual_evaluate.R


job_file="/gpfs/gibbs/pi/zhao/lx94/SWIFT/code/simulation/residual/EUR_residual_analysis.txt"
> $job_file  # Empty the job file if it already exists

pop=EUR
sample_size=ukbb
i=1

for sim_i in {1..5}; do
for repeat_i in {1..4}; do
for p in 0.1 0.01 0.001 5e-04; do

out_check="/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/residual/${pop}_sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${sample_size}_prune_snplist_${i}_${pop}_tune_GWAS_residual_approxTRUE_ratio3.00_repeat${repeat_i}.txt"

if [[ ! -e "${out_check}" ]]; then

echo "module load miniconda; conda activate r_env; Rscript --vanilla /gpfs/gibbs/pi/zhao/lx94/SWIFT/code/simulation/residual/EUR_residual_evaluate.R ${sim_i} ${repeat_i} ${p}" >> $job_file

fi

done
done
done


module load dSQ
cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/code/simulation/residual/

dsq --job-file /gpfs/gibbs/pi/zhao/lx94/SWIFT/code/simulation/residual/EUR_residual_analysis.txt --partition=scavenge,day --requeue --mem=30G --cpus-per-task=1 --ntasks=1 --nodes=1 --time=24:00:00 --mail-type=ALL
sbatch dsq-EUR_residual_analysis-$(date +%Y-%m-%d).sh
