## Step1: Format both real and subsample GWAS into method format
## real GWAS
library(data.table)

for (trait in c("WBC","NEU","LYM","MON","EOS","RBC","HCT","MCH","MCV","HB","ALT","ALP","GGT")){

pop = "EUR"

## full snplist
sumstat_data_full_snplist = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/summary_data/clean/",trait, "_", pop, "_inter_clean_full_snplist.txt"))

PRScsx_data_full_snplist = sumstat_data_full_snplist[,c("SNP","A1","A2","BETA","P")]
write.table(PRScsx_data_full_snplist, file=paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/summary_data/PRScsx/",trait, "_", pop, "_inter_PRScsx_full_snplist.txt"), 
            row.names=F, col.names=T, quote=F, append=F, sep = "\t")

sumstat_SDPRX_full_snplist = sumstat_data_full_snplist[,c("SNP","A1","A2","Z","P","N")]
write.table(sumstat_SDPRX_full_snplist, file=paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/summary_data/SDPRX/",trait,"_",pop,"_inter_SDPRX_full_snplist.txt"), 
            row.names=F, col.names=T, quote=F, append=F, sep = "\t")


## prune snplist
for (i in c(1:4)){

pop2 = "EAS"

sumstat_data_prune_snplist = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/summary_data/clean/",trait, "_", pop, "_inter_clean_",pop2,"_prune_snplist_",i,".txt"))

PRScsx_data_prune_snplist = sumstat_data_prune_snplist[,c("SNP","A1","A2","BETA","P")]
write.table(PRScsx_data_prune_snplist, file=paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/summary_data/PRScsx/",trait, "_", pop, "_inter_PRScsx_",pop2,"_prune_snplist_",i,".txt"), 
            row.names=F, col.names=T, quote=F, append=F, sep = "\t")

sumstat_SDPRX_prune_snplist = sumstat_data_prune_snplist[,c("SNP","A1","A2","Z","P","N")]
write.table(sumstat_SDPRX_prune_snplist, file=paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/summary_data/SDPRX/",trait,"_",pop,"_inter_SDPRX_",pop2,"_prune_snplist_",i,".txt"), 
            row.names=F, col.names=T, quote=F, append=F, sep = "\t")

}
}


## subsample GWAS
## train
library(data.table)

for (trait in c("WBC","NEU","LYM","MON","EOS","RBC","HCT","MCH","MCV","HB","ALT","ALP","GGT")){

pop = "EAS"

for (rpt in c(1:4)){
## full snplist
sumstat_data_full_snplist = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/summary_data/subsample/clean/",trait,"_full_snplist_",pop,"_train_GWAS_approxFALSE_ratio3.00_repeat",rpt,".txt"))

PRScsx_data_full_snplist = sumstat_data_full_snplist[,c("SNP","A1","A2","BETA","P")]
write.table(PRScsx_data_full_snplist, file=paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/summary_data/subsample/PRScsx/",trait,"_full_snplist_",pop,"_train_PRScsx_approxFALSE_ratio3.00_repeat",rpt,".txt"), 
            row.names=F, col.names=T, quote=F, append=F, sep = "\t")

sumstat_SDPRX_full_snplist = sumstat_data_full_snplist[,c("SNP","A1","A2","Z","P","N")]
write.table(sumstat_SDPRX_full_snplist, file=paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/summary_data/subsample/SDPRX/",trait,"_full_snplist_",pop,"_train_SDPRX_approxFALSE_ratio3.00_repeat",rpt,".txt"), 
            row.names=F, col.names=T, quote=F, append=F, sep = "\t")


## prune snplist
for (i in c(1:4)){
for (approx in c("TRUE","FALSE")){

sumstat_data_prune_snplist = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/summary_data/subsample/clean/",trait,"_prune_snplist_",i,"_",pop,"_train_GWAS_approx",approx,"_ratio3.00_repeat",rpt,".txt"))

PRScsx_data_prune_snplist = sumstat_data_prune_snplist[,c("SNP","A1","A2","BETA","P")]
write.table(PRScsx_data_prune_snplist, file=paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/summary_data/subsample/PRScsx/",trait,"_prune_snplist_",i,"_",pop,"_train_PRScsx_approx",approx,"_ratio3.00_repeat",rpt,".txt"), 
            row.names=F, col.names=T, quote=F, append=F, sep = "\t")

sumstat_SDPRX_prune_snplist = sumstat_data_prune_snplist[,c("SNP","A1","A2","Z","P","N")]
write.table(sumstat_SDPRX_prune_snplist, file=paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/summary_data/subsample/SDPRX/",trait,"_prune_snplist_",i,"_",pop,"_train_SDPRX_approx",approx,"_ratio3.00_repeat",rpt,".txt"), 
            row.names=F, col.names=T, quote=F, append=F, sep = "\t")

}
}

}
}

## tune
library(data.table)

for (trait in c("WBC","NEU","LYM","MON","EOS","RBC","HCT","MCH","MCV","HB","ALT","ALP","GGT")){

pop = "EAS"

for (rpt in c(1:4)){
## full snplist
sumstat_data_full_snplist = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/summary_data/subsample/clean/",trait,"_full_snplist_",pop,"_tune_GWAS_approxFALSE_ratio3.00_repeat",rpt,".txt"))

sumstat_MIX_full_snplist = sumstat_data_full_snplist[,c("SNP","A1","A2","BETA","SE","Z","P","N")]
write.table(sumstat_MIX_full_snplist, file=paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/summary_data/subsample/MIX/",trait,"_full_snplist_",pop,"_tune_MIX_approxFALSE_ratio3.00_repeat",rpt,".txt"), 
            row.names=F, col.names=T, quote=F, append=F, sep = "\t")


## prune snplist
for (i in c(1:4)){
for (approx in c("TRUE","FALSE")){

sumstat_data_prune_snplist = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/summary_data/subsample/clean/",trait,"_prune_snplist_",i,"_",pop,"_tune_GWAS_approx",approx,"_ratio3.00_repeat",rpt,".txt"))

sumstat_MIX_prune_snplist = sumstat_data_prune_snplist[,c("SNP","A1","A2","BETA","SE","Z","P","N")]
write.table(sumstat_MIX_prune_snplist, file=paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/summary_data/subsample/MIX/",trait,"_prune_snplist_",i,"_",pop,"_tune_MIX_approx",approx,"_ratio3.00_repeat",rpt,".txt"), 
            row.names=F, col.names=T, quote=F, append=F, sep = "\t")

}
}

}
}