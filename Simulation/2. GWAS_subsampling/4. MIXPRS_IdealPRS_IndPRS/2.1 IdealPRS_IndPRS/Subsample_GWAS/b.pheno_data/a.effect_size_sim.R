# 0. Inter all pop snp list
library(data.table)

EUR_snplist = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/snplist_data/EUR_inter_snplist.txt"),header=F)
EAS_snplist = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/snplist_data/EAS_inter_snplist.txt"),header=F)
AFR_snplist = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/snplist_data/AFR_inter_snplist.txt"),header=F)
SAS_snplist = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/snplist_data/SAS_inter_snplist.txt"),header=F)
AMR_snplist = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/snplist_data/AMR_inter_snplist.txt"),header=F)

snplist = intersect(EUR_snplist$V1,EAS_snplist$V1)
snplist = intersect(snplist,AFR_snplist$V1)
snplist = intersect(snplist,SAS_snplist$V1)
snplist = intersect(snplist,AMR_snplist$V1)

snp = read.table("/gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/info/snp_rs_id_infor_hm3.snplist")
snp = snp[which(snp$V1 %in% snplist),]
snp = data.table(SNP = snp)

write.table(snp, 
            file="/gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/info/snp_allpop_inter_hm3.snplist", 
            append = F, quote = F, sep = "\t", row.names = F,
            col.names = T)

## 1. True effect size simulation for all SNPs
library(MASS) # mvrnorm
library(data.table)
library(tidyr)
library(dplyr)

# set up
n_pop = 5
rep_num = 5
h2 = rep(0.4, n_pop) # heritability
h2_sqrt = sqrt(h2)
rhog = 0.8

pop_list = c("EUR","EAS","AFR","SAS","AMR")
p_list = c(0.001,0.01,0.0005,0.1)

save_path <- paste0("/gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/effect_data/")
if (!dir.exists(save_path)) {
  dir.create(save_path, recursive = TRUE)
}

# obtain common snp
snp_info = fread("/gpfs/gibbs/pi/zhao/xz674/data/sim_data/geno_data/info/snp_infor_mega+hm3")
snp_info = snp_info[complete.cases(snp_info),]
snp = read.table("/gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/info/snp_rs_id_infor_hm3.snplist")
snp_info = snp_info[which(snp_info$rs_id %in% snp$V1),]

snp_split <- snp_info %>% separate(SNP,into=c("rsid","position2","noncoding","coding"),sep=":")
snp = snp_info$rs_id
S = length(snp)

write.table(snp, 
            file="/gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/info/snp_common_hm3.snplist", 
            append = F, quote = F, sep = "\t", row.names = F,
            col.names = F)

EUR_snp_infor = snp_split[,c("rsid","noncoding","coding","FREQ_A1_EUR")]
EAS_snp_infor = snp_split[,c("rsid","noncoding","coding","FREQ_A1_EAS")]
AFR_snp_infor = snp_split[,c("rsid","noncoding","coding","FREQ_A1_AFR")]
SAS_snp_infor = snp_split[,c("rsid","noncoding","coding","FREQ_A1_SAS")]
AMR_snp_infor = snp_split[,c("rsid","noncoding","coding","FREQ_A1_AMR")]

colnames(EUR_snp_infor) = c("rsid","noncoding","coding","FREQ_coding")
colnames(EAS_snp_infor) = c("rsid","noncoding","coding","FREQ_coding")
colnames(AFR_snp_infor) = c("rsid","noncoding","coding","FREQ_coding")
colnames(SAS_snp_infor) = c("rsid","noncoding","coding","FREQ_coding")
colnames(AMR_snp_infor) = c("rsid","noncoding","coding","FREQ_coding")

write.table(EUR_snp_infor, 
            file="/gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/info/EUR_snp_infor", 
            append = F, quote = F, sep = "\t", row.names = F,
            col.names = T)
write.table(EAS_snp_infor, 
            file="/gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/info/EAS_snp_infor", 
            append = F, quote = F, sep = "\t", row.names = F,
            col.names = T)
write.table(AFR_snp_infor, 
            file="/gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/info/AFR_snp_infor", 
            append = F, quote = F, sep = "\t", row.names = F,
            col.names = T)
write.table(SAS_snp_infor, 
            file="/gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/info/SAS_snp_infor", 
            append = F, quote = F, sep = "\t", row.names = F,
            col.names = T)
write.table(AMR_snp_infor, 
            file="/gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/info/AMR_snp_infor", 
            append = F, quote = F, sep = "\t", row.names = F,
            col.names = T)

# find common SNPs and assign effect SNPs among them
snp_common = snp_info[which(snp_info$FREQ_A1_EUR >= 0.01 & snp_info$FREQ_A1_EUR <= 0.99),]
snp_common = snp_common[which(snp_common$FREQ_A1_EAS >= 0.01 & snp_common$FREQ_A1_EAS <= 0.99),]
snp_common = snp_common[which(snp_common$FREQ_A1_AFR >= 0.01 & snp_common$FREQ_A1_AFR <= 0.99),]
snp_common = snp_common[which(snp_common$FREQ_A1_SAS >= 0.01 & snp_common$FREQ_A1_SAS <= 0.99),]
snp_common = snp_common[which(snp_common$FREQ_A1_AMR >= 0.01 & snp_common$FREQ_A1_AMR <= 0.99),]

# intersect with 1kg (1009128 -> 905640)
for (pop in c("EUR","EAS","AFR","SAS","AMR")){
  if (pop == "EUR"){ssize = "100K"}else{ssize = "20K"}
  for (geno in c(ssize, "1kg")){
    bim_data = read.table(paste0("/gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/",pop,"/discover_validate/",pop,"_",geno,"_all.bim"))
    colnames(bim_data) = c("chr","snpid","nonthing","position","minor","major")
    snp_common = snp_common[which(snp_common$rs_id %in% bim_data$snpid),]
  }
}

# intersect with ukbb (899505)
bim_data = read.table("/gpfs/gibbs/pi/zhao/xz674/data/ukbb_data/geno_data/toy_example.bim")
colnames(bim_data) = c("chr","snpid","nonthing","position","minor","major")
snp_common = snp_common[which(snp_common$rs_id %in% bim_data$snpid),]

write.table(snp_common$rs_id, 
            file="/gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/info/snp_common_allpop_inter_hm3.snplist", 
            append = F, quote = F, sep = "\t", row.names = F,
            col.names = F)

print("Start!")
set.seed(0)
for (p in p_list){
  # need to obtain the minor allele effect
  num_select = round(p * S)
  p_variable <- rep(0, S)
  selected_indices <- sample(which(snp %in% snp_common$rs_id), size = num_select, replace = FALSE)
  p_variable[selected_indices] <- 1
  
  rhog_matrix <- diag(x = 1 - rhog, nrow = n_pop) + matrix(rep(rhog, n_pop ^ 2), nrow = n_pop)
  covg_matrix <- 1 / (p * S) * diag(h2_sqrt) %*% rhog_matrix %*% diag(h2_sqrt)
  
  for (i in 1:rep_num){
    print(paste0("Start: p=", p, ", sim=", i))
    
    while (TRUE) {
      beta_true <- data.frame(rep(0, S))
      for (pop in pop_list){
        beta_true[[paste0(pop, "_beta_true")]] = rep(0, S)
      }
      beta_true <- beta_true[, -1]
      beta_true[,colnames(beta_true)] <- mvrnorm(S, rep(0, n_pop), Sigma = covg_matrix)
      for (pop in pop_list){
        beta_true[,paste0(pop, "_beta_true")] = p_variable * beta_true[,paste0(pop, "_beta_true")]
      }
      
      if (all(abs(colSums(beta_true ^ 2) - h2) < 1e-2)) {break}
    }
    
    for (pop in pop_list){
      print(paste0(pop, " h2: ",sum(beta_true[[paste0(pop, "_beta_true")]] ^ 2)))
    }
    
    print(paste0("EUR_EAS corr: ",cor(beta_true$EUR_beta_true,beta_true$EAS_beta_true)))
    print(paste0("EUR_AFR corr: ",cor(beta_true$EUR_beta_true,beta_true$AFR_beta_true)))
    print(paste0("EUR_SAS corr: ",cor(beta_true$EUR_beta_true,beta_true$SAS_beta_true)))
    print(paste0("EUR_AMR corr: ",cor(beta_true$EUR_beta_true,beta_true$AMR_beta_true)))
    
    for (pop in pop_list){
      trait_causal_origin = data.table(snp_id = snp,
                                       coding = snp_split$coding,
                                       effect = beta_true[[paste0(pop, "_beta_true")]])

      # simulated geno allele benchmark
      bim_name = paste0("/gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/",pop,"/discover_validate/",pop,"_100K_all.bim")
      trait_bim.infor <- fread(bim_name)
      colnames(trait_bim.infor) <- c("chr","snpid","nonthing","position","minor","major")
      trait_bim.infor = trait_bim.infor[,c("snpid","minor")]
      trait_causal = trait_causal_origin[trait_bim.infor, on = .(snp_id = snpid)]
      trait_causal = na.omit(trait_causal)
      trait_causal$effect <- ifelse(trait_causal$coding!=trait_causal$minor,-trait_causal$effect,trait_causal$effect)
      trait_causal = trait_causal[which(trait_causal$effect != 0),]
      # write effect size
      write.table(trait_causal[,c(1,3)], 
                  file=paste0(save_path, pop,"_sim",i,"_h2",h2[1],"_p",p,"_rhog",rhog,"_simu.txt"), 
                  append = F, quote = F, sep = "\t", row.names = F, col.names = F)
      print(paste0("Write ",pop," effect size for simulated geno."))
      # write causal snps
      write.table(trait_causal[,c(1)], 
                  file=paste0(save_path, pop,"_sim",i,"_h2",h2[1],"_p",p,"_rhog",rhog,"_simu.snplist"), 
                  append = F, quote = F, sep = "\t", row.names = F, col.names = F)
      print(paste0("Write ",pop," snplist for simulated geno."))

      # 1kg allele benchmark
      bim_name = paste0("/gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/",pop,"/discover_validate/",pop,"_1kg_all.bim")
      trait_bim.infor <- fread(bim_name)
      colnames(trait_bim.infor) <- c("chr","snpid","nonthing","position","minor","major")
      trait_bim.infor = trait_bim.infor[,c("snpid","minor")]
      trait_causal = trait_causal_origin[trait_bim.infor, on = .(snp_id = snpid)]
      trait_causal = na.omit(trait_causal)
      trait_causal$effect <- ifelse(trait_causal$coding!=trait_causal$minor,-trait_causal$effect,trait_causal$effect)
      trait_causal = trait_causal[which(trait_causal$effect != 0),]
      # write effect size
      write.table(trait_causal[,c(1,3)], 
                  file=paste0(save_path, pop,"_sim",i,"_h2",h2[1],"_p",p,"_rhog",rhog,"_1kg.txt"), 
                  append = F, quote = F, sep = "\t", row.names = F, col.names = F)
      print(paste0("Write ",pop," effect size for 1kg geno."))
      # write causal snps
      write.table(trait_causal[,c(1)], 
                  file=paste0(save_path, pop,"_sim",i,"_h2",h2[1],"_p",p,"_rhog",rhog,"_1kg.snplist"), 
                  append = F, quote = F, sep = "\t", row.names = F, col.names = F)
      print(paste0("Write ",pop," snplist for 1kg geno."))

      # ukbb allele benchmark
      if (pop == "EUR"){
        bim_name = paste0("/gpfs/gibbs/pi/zhao/xz674/data/ukbb_data/geno_data/toy_example.bim")
        trait_bim.infor <- fread(bim_name)
        colnames(trait_bim.infor) <- c("chr","snpid","nonthing","position","minor","major")
        trait_bim.infor = trait_bim.infor[,c("snpid","minor")]
        trait_causal = trait_causal_origin[trait_bim.infor, on = .(snp_id = snpid)]
        trait_causal = na.omit(trait_causal)
        trait_causal$effect <- ifelse(trait_causal$coding!=trait_causal$minor,-trait_causal$effect,trait_causal$effect)
        trait_causal = trait_causal[which(trait_causal$effect != 0),]
        # write effect size
        write.table(trait_causal[,c(1,3)], 
                    file=paste0(save_path, pop,"_sim",i,"_h2",h2[1],"_p",p,"_rhog",rhog,"_ukbb.txt"), 
                    append = F, quote = F, sep = "\t", row.names = F, col.names = F)
        print(paste0("Write ",pop," effect size for ukbb geno."))
        # write causal snps
        write.table(trait_causal[,c(1)], 
                    file=paste0(save_path, pop,"_sim",i,"_h2",h2[1],"_p",p,"_rhog",rhog,"_ukbb.snplist"), 
                    append = F, quote = F, sep = "\t", row.names = F, col.names = F)
        print(paste0("Write ",pop," snplist for ukbb geno."))
      }
    }
  }
}

print("Finish!")
