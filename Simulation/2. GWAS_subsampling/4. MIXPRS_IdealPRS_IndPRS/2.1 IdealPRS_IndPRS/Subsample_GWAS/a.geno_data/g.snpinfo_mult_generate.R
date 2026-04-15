library(data.table)
library(dplyr)

# 5 populations
pop_list <- c("AFR", "AMR", "EAS", "EUR", "SAS")

# generate mult snpinfo function
generate_snpinfo_mult <- function(geno, geno_type, split, cv){
  # Populations data frames
  snpinfo_mult <- list()
  for (pop in pop_list){
    file_name <- paste0("/gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/", pop, "/", split, "/", pop, "_", geno, "_", cv, ".bim")
    snpinfo_mult[[pop]] <- read.table(file_name, header = F)[,-3]
  }
  df <- full_join(snpinfo_mult[[pop_list[1]]], snpinfo_mult[[pop_list[2]]], by=c("V1","V2","V4","V5","V6"))
  for (i in 3:length(pop_list)){
    df <- df %>% full_join(snpinfo_mult[[pop_list[i]]],by=c("V1","V2","V4","V5","V6"))
  }
  
  # Merge and calculate frequencies and flip status
  colnames(df) <- c("CHR","SNP","BP","A1","A2")
  
  for (pop in pop_list){
    if (split == "discover_validate"){
      file_name <- paste0("/gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/ref_data/", split, "/", geno, "/ldblk_", geno_type, "_", tolower(pop), "/snpinfo_", geno_type, "_hm3")
    } else{
      file_name <- paste0("/gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/ref_data/", split, "/", geno, "_", cv, "/ldblk_", geno_type, "_", tolower(pop), "/snpinfo_", geno_type, "_hm3")
    }
    maf_pop <- read.table(file_name, header = T)
    df <- full_join(df, maf_pop, by = c("CHR", "SNP", "BP"), suffix = c("", "_pop"))
    df <- df %>% 
      mutate(
        FRQ = case_when(
          A1 == A1_pop & A2 == A2_pop ~ MAF ,
          A1 == A2_pop & A2 == A1_pop ~ 1 - MAF ,
          TRUE ~ NA_real_
        )
      )
    names(df)[names(df) == "FRQ"] <- paste0("FRQ_", pop)
    df <- df %>%
      select(-A1_pop, -A2_pop, -MAF)
  }
  
  df[is.na(df)] = 0
  p <- nrow(df)
  
  for (pop in pop_list){
    df <- cbind(df, setNames(data.frame(rep(1, p)), paste0("FLP_", pop)))
  }
  df <- df[order(df$CHR,df$BP),]
  
  if (split == "discover_validate"){
    file_name <- paste0("/gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/ref_data/", split, "/", geno, "/snpinfo_mult_", geno_type, "_hm3")
  } else{
    file_name <- paste0("/gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/ref_data/", split, "/", geno, "_", cv, "/snpinfo_mult_", geno_type, "_hm3")
  }
  
  print(paste0("write multi snpinfo for ", geno, "_", split, "_", cv))
  write.table(df, file_name, sep = "\t", row.names = FALSE, quote = F)
}

# main
for (geno in c("100K")){ # "1kg", "1kg1", "1kg2"
  geno_type <- ifelse(geno %in% c("1kg", "1kg1", "1kg2"), "1kg", geno)
  for (split in c("discover_validate")){ # "discover", "validate", 
    if (split == "discover_validate"){
      cv <- "all"
      generate_snpinfo_mult(geno, geno_type, split, cv)
    } else{
      for (ff in 1:4){
        cv <- paste0("fold", ff)
        generate_snpinfo_mult(geno, geno_type, split, cv)
      }
    }
  }
}

