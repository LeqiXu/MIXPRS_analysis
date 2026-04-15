## 1. Generated Data For Each Method
library(data.table)
library(stringr)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

pop = args[1]

update_name = fread("/gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/info/snp_map_hm3.snplist", header=F)
snp_list = fread("/gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/info/snp_allpop_inter_hm3.snplist")

# set up
n_fold = 4
h2 = 0.4
rhog = 0.8

print("Start!")

snp_info = fread(paste0("/gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/info/",pop,"_snp_infor"))
# if (pop == "EUR"){ssize = "100K"}else{ssize = c("20K", "100K")}
if (pop == "EUR"){ssize = "ukbb"}else{ssize = c("20K", "100K")}
# geno_list = c(ssize, "1kg", "1kg1", "1kg2")
geno_list = c(ssize)
for (p in c(0.001,0.01,5e-04,0.1)){
  for (i in 1:5){
    print(paste0("pop=", pop, ", i=", i, ", p=", p, ", rhog=", rhog))
    
    # ## discover_validate
    # for (geno in geno_list){
    #   if (! file.exists(paste0("/gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/summary_data/",pop,"/discover_validate/clean/",pop,"_sim",i,"_h2",h2,"_p",p,"_rhog",rhog,"_",geno,"_clean_real_all.txt"))){
    #     print(paste0("split=discover_validate, geno=", geno))
    #     data = fread(paste0("/gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/summary_data/",pop,"/discover_validate/",pop,"_sim",i,"_h2",h2,"_p",p,"_rhog",rhog,"_",geno,"_all.assoc.linear"))
    #     data = data[which(data$SNP %in% snp_list$SNP),]
    #     data = data[snp_info, on = .(SNP = rsid)]
    #     data$A2 = ifelse(data$A1 == snp_info$noncoding, snp_info$coding,snp_info$noncoding)
    #     data$MAF = ifelse(data$A1 == snp_info$noncoding, 1 - snp_info$FREQ_coding, snp_info$FREQ_coding)
    #     data = data[which(data$SNP %in% snp_list$SNP),]
        
    #     # clean
    #     clean_data = data[,c("SNP","CHR","BP","A1","A2","NMISS","MAF","BETA","STAT","P")]  # NMISS? NCHROB? 
    #     colnames(clean_data) = c("SNP","CHR","POS","A1","A2","N","MAF","BETA","Z","P")
    #     clean_data$SE = clean_data$BETA / clean_data$Z
    #     clean_data = clean_data[update_name, on = .(SNP = V2)]
    #     colnames(clean_data) = c("SNP","CHR","POS","A1","A2","N","MAF","BETA","Z","P","SE","rsid")
    #     clean_data = na.omit(clean_data)
    #     write.table(clean_data, file=paste0("/gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/summary_data/",pop,"/discover_validate/clean/",pop,"_sim",i,"_h2",h2,"_p",p,"_rhog",rhog,"_",geno,"_clean_real_all.txt"), 
    #                 row.names=F, col.names=T, quote=F, append=F, sep = "\t") 
    #   }
    # }
    
    for (ff in 1:n_fold){
      ## discover
      for (geno in geno_list){
        if (! file.exists(paste0("/gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/summary_data/clean/",pop,"_sim",i,"_h2",h2,"_p",p,"_rhog",rhog,"_",geno,"_clean_real_fold",ff,".txt"))){
          print(paste0("split=discover, geno=", geno, ", fold=", ff))
          data = fread(paste0("/gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/summary_data/",pop,"/discover/",pop,"_sim",i,"_h2",h2,"_p",p,"_rhog",rhog,"_",geno,"_fold",ff,".assoc.linear"))
          data = data[which(data$SNP %in% snp_list$SNP),]
          data = data[snp_info, on = .(SNP = rsid)]
          data$A2 = ifelse(data$A1 == snp_info$noncoding, snp_info$coding,snp_info$noncoding)
          data$MAF = ifelse(data$A1 == snp_info$noncoding, 1 - snp_info$FREQ_coding, snp_info$FREQ_coding)
          data = data[which(data$SNP %in% snp_list$SNP),]
          
          # clean
          clean_data = data[,c("SNP","CHR","BP","A1","A2","NMISS","MAF","BETA","STAT","P")]
          colnames(clean_data) = c("SNP","CHR","POS","A1","A2","N","MAF","BETA","Z","P")
          clean_data$SE = clean_data$BETA / clean_data$Z
          clean_data = clean_data[update_name, on = .(SNP = V2)]
          colnames(clean_data) = c("SNP","CHR","POS","A1","A2","N","MAF","BETA","Z","P","SE","rsid")
          clean_data = na.omit(clean_data)
          write.table(clean_data, file=paste0("/gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/summary_data/clean/",pop,"_sim",i,"_h2",h2,"_p",p,"_rhog",rhog,"_",geno,"_clean_real_fold",ff,".txt"), 
                      row.names=F, col.names=T, quote=F, append=F, sep = "\t") 
        }
      }
      
      # ## validate
      # for (geno in geno_list){
      #   if (! file.exists(paste0("/gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/summary_data/",pop,"/validate/clean/",pop,"_sim",i,"_h2",h2,"_p",p,"_rhog",rhog,"_",geno,"_clean_real_fold",ff,".txt"))){
      #     print(paste0("split=validate, geno=", geno, ", fold=", ff))
      #     data = fread(paste0("/gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/summary_data/",pop,"/validate/",pop,"_sim",i,"_h2",h2,"_p",p,"_rhog",rhog,"_",geno,"_fold",ff,".assoc.linear"))
      #     data = data[which(data$SNP %in% snp_list$SNP),]
      #     data = data[snp_info, on = .(SNP = rsid)]
      #     data$A2 = ifelse(data$A1 == snp_info$noncoding, snp_info$coding,snp_info$noncoding)
      #     data$MAF = ifelse(data$A1 == snp_info$noncoding, 1 - snp_info$FREQ_coding, snp_info$FREQ_coding)
      #     data = data[which(data$SNP %in% snp_list$SNP),]
          
      #     # clean
      #     clean_data = data[,c("SNP","CHR","BP","A1","A2","NMISS","MAF","BETA","STAT","P")]
      #     colnames(clean_data) = c("SNP","CHR","POS","A1","A2","N","MAF","BETA","Z","P")
      #     clean_data$SE = clean_data$BETA / clean_data$Z
      #     clean_data = clean_data[update_name, on = .(SNP = V2)]
      #     colnames(clean_data) = c("SNP","CHR","POS","A1","A2","N","MAF","BETA","Z","P","SE","rsid")
      #     clean_data = na.omit(clean_data)
      #     write.table(clean_data, file=paste0("/gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/summary_data/",pop,"/validate/clean/",pop,"_sim",i,"_h2",h2,"_p",p,"_rhog",rhog,"_",geno,"_clean_real_fold",ff,".txt"), 
      #                 row.names=F, col.names=T, quote=F, append=F, sep = "\t") 
      #   }
      # }
    }
  }
}

# # UKBB data
# geno = "ukbb"
# if (pop == "EUR"){
#   snp_info = fread(paste0("/gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/info/",pop,"_snp_infor"))
#   for (p in c(0.001,0.01,5e-04,0.1)){
#     for (i in 1:5){
#       print(paste0("pop=", pop, ", i=", i, ", p=", p, ", rhog=", rhog))
      
#       # ## discover_validate
#       # if (! file.exists(paste0("/gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/summary_data/",pop,"/discover_validate/clean/",pop,"_sim",i,"_h2",h2,"_p",p,"_rhog",rhog,"_",geno,"_clean_real_all.txt"))){
#       #   print(paste0("split=discover_validate, geno=", geno))
#       #   data = fread(paste0("/gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/summary_data/",pop,"/discover_validate/",pop,"_sim",i,"_h2",h2,"_p",p,"_rhog",rhog,"_",geno,"_all.assoc.linear"))
#       #   data = data[which(data$SNP %in% snp_list$SNP),]
#       #   data = data[snp_info, on = .(SNP = rsid)]
#       #   data$A2 = ifelse(data$A1 == snp_info$noncoding, snp_info$coding,snp_info$noncoding)
#       #   data$MAF = ifelse(data$A1 == snp_info$noncoding, 1 - snp_info$FREQ_coding, snp_info$FREQ_coding)
#       #   data = data[which(data$SNP %in% snp_list$SNP),]

#       #   # clean
#       #   clean_data = data[,c("SNP","CHR","BP","A1","A2","NMISS","MAF","BETA","STAT","P")]  # NMISS? NCHROB?
#       #   colnames(clean_data) = c("SNP","CHR","POS","A1","A2","N","MAF","BETA","Z","P")
#       #   clean_data$SE = clean_data$BETA / clean_data$Z
#       #   clean_data = clean_data[update_name, on = .(SNP = V2)]
#       #   colnames(clean_data) = c("SNP","CHR","POS","A1","A2","N","MAF","BETA","Z","P","SE","rsid")
#       #   clean_data = na.omit(clean_data)
#       #   write.table(clean_data, file=paste0("/gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/summary_data/",pop,"/discover_validate/clean/",pop,"_sim",i,"_h2",h2,"_p",p,"_rhog",rhog,"_",geno,"_clean_real_all.txt"),
#       #               row.names=F, col.names=T, quote=F, append=F, sep = "\t")
#       # }
      
#       ## discover
#       if (! file.exists(paste0("/gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/summary_data/clean/",pop,"_sim",i,"_h2",h2,"_p",p,"_rhog",rhog,"_",geno,"_clean_real.txt"))){
#         print(paste0("split=discover, geno=", geno))
#         data = fread(paste0("/gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/summary_data/",pop,"/discover/",pop,"_sim",i,"_h2",h2,"_p",p,"_rhog",rhog,"_",geno,".assoc.linear"))
#         data = data[which(data$SNP %in% snp_list$SNP),]
#         data = data[snp_info, on = .(SNP = rsid)]
#         data$A2 = ifelse(data$A1 == snp_info$noncoding, snp_info$coding,snp_info$noncoding)
#         data$MAF = ifelse(data$A1 == snp_info$noncoding, 1 - snp_info$FREQ_coding, snp_info$FREQ_coding)
#         data = data[which(data$SNP %in% snp_list$SNP),]
        
#         # clean
#         clean_data = data[,c("SNP","CHR","BP","A1","A2","NMISS","MAF","BETA","STAT","P")]
#         colnames(clean_data) = c("SNP","CHR","POS","A1","A2","N","MAF","BETA","Z","P")
#         clean_data$SE = clean_data$BETA / clean_data$Z
#         clean_data = clean_data[update_name, on = .(SNP = V2)]
#         colnames(clean_data) = c("SNP","CHR","POS","A1","A2","N","MAF","BETA","Z","P","SE","rsid")
#         clean_data = na.omit(clean_data)
#         write.table(clean_data, file=paste0("/gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/summary_data/clean/",pop,"_sim",i,"_h2",h2,"_p",p,"_rhog",rhog,"_",geno,"_clean_real.txt"), 
#                     row.names=F, col.names=T, quote=F, append=F, sep = "\t") 
#       }

      # ## validate
      # if (! file.exists(paste0("/gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/summary_data/",pop,"/validate/clean/",pop,"_sim",i,"_h2",h2,"_p",p,"_rhog",rhog,"_",geno,"_clean_real.txt"))){
      #   print(paste0("split=validate, geno=", geno))
      #   data = fread(paste0("/gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/summary_data/",pop,"/validate/",pop,"_sim",i,"_h2",h2,"_p",p,"_rhog",rhog,"_",geno,".assoc.linear"))
      #   data = data[which(data$SNP %in% snp_list$SNP),]
      #   data = data[snp_info, on = .(SNP = rsid)]
      #   data$A2 = ifelse(data$A1 == snp_info$noncoding, snp_info$coding,snp_info$noncoding)
      #   data$MAF = ifelse(data$A1 == snp_info$noncoding, 1 - snp_info$FREQ_coding, snp_info$FREQ_coding)
      #   data = data[which(data$SNP %in% snp_list$SNP),]
        
      #   # clean
      #   clean_data = data[,c("SNP","CHR","BP","A1","A2","NMISS","MAF","BETA","STAT","P")]
      #   colnames(clean_data) = c("SNP","CHR","POS","A1","A2","N","MAF","BETA","Z","P")
      #   clean_data$SE = clean_data$BETA / clean_data$Z
      #   clean_data = clean_data[update_name, on = .(SNP = V2)]
      #   colnames(clean_data) = c("SNP","CHR","POS","A1","A2","N","MAF","BETA","Z","P","SE","rsid")
      #   clean_data = na.omit(clean_data)
      #   write.table(clean_data, file=paste0("/gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/summary_data/",pop,"/validate/clean/",pop,"_sim",i,"_h2",h2,"_p",p,"_rhog",rhog,"_",geno,"_clean_real.txt"), 
      #               row.names=F, col.names=T, quote=F, append=F, sep = "\t") 
      # }
    }
  }
}

print("Finish!")