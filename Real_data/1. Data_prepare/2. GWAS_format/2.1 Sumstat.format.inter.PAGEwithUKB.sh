## Step1: Format both real and subsample GWAS into method format
## real GWAS
library(data.table)

for (trait in c("Height","BMI","SBP","DBP","PLT")){

for (pop in c("AFR")){

## prune snplist
i=1
for (pop2 in c("AFR")){

if (pop2 != pop){
sumstat_data_prune_snplist = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/summary_data/clean/",trait, "_", pop, "_inter_UKB_clean_",pop2,"_prune_snplist_",i,".txt"))

PRScsx_data_prune_snplist = sumstat_data_prune_snplist[,c("SNP","A1","A2","BETA","P")]
write.table(PRScsx_data_prune_snplist, file=paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/summary_data/PRScsx/",trait, "_", pop, "_inter_UKB_PRScsx_",pop2,"_prune_snplist_",i,".txt"), 
            row.names=F, col.names=T, quote=F, append=F, sep = "\t")

sumstat_SDPRX_prune_snplist = sumstat_data_prune_snplist[,c("SNP","A1","A2","Z","P","N")]
write.table(sumstat_SDPRX_prune_snplist, file=paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/summary_data/SDPRX/",trait,"_",pop,"_inter_UKB_SDPRX_",pop2,"_prune_snplist_",i,".txt"), 
            row.names=F, col.names=T, quote=F, append=F, sep = "\t")

}
}

}

}


## subsample GWAS
## train
library(data.table)

for (trait in c("Height","BMI","SBP","DBP","PLT")){

for (pop in c("AFR")){

for (rpt in c(1:4)){

## prune snplist
i=1
approx="TRUE"


sumstat_data_prune_snplist = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/summary_data/subsample/clean/",trait,"_UKB_prune_snplist_",i,"_",pop,"_train_GWAS_approx",approx,"_ratio3.00_repeat",rpt,".txt"))

PRScsx_data_prune_snplist = sumstat_data_prune_snplist[,c("SNP","A1","A2","BETA","P")]
write.table(PRScsx_data_prune_snplist, file=paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/summary_data/subsample/PRScsx/",trait,"_UKB_prune_snplist_",i,"_",pop,"_train_PRScsx_approx",approx,"_ratio3.00_repeat",rpt,".txt"), 
            row.names=F, col.names=T, quote=F, append=F, sep = "\t")

sumstat_SDPRX_prune_snplist = sumstat_data_prune_snplist[,c("SNP","A1","A2","Z","P","N")]
write.table(sumstat_SDPRX_prune_snplist, file=paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/summary_data/subsample/SDPRX/",trait,"_UKB_prune_snplist_",i,"_",pop,"_train_SDPRX_approx",approx,"_ratio3.00_repeat",rpt,".txt"), 
            row.names=F, col.names=T, quote=F, append=F, sep = "\t")

}
}

}

## tune
library(data.table)

for (trait in c("Height","BMI","SBP","DBP","PLT")){

for (pop in c("AFR")){

for (rpt in c(1:4)){

## prune snplist
i=1
approx="TRUE"

sumstat_data_prune_snplist = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/summary_data/subsample/clean/",trait,"_UKB_prune_snplist_",i,"_",pop,"_tune_GWAS_approx",approx,"_ratio3.00_repeat",rpt,".txt"))

sumstat_MIX_prune_snplist = sumstat_data_prune_snplist[,c("SNP","A1","A2","BETA","SE","Z","P","N")]
write.table(sumstat_MIX_prune_snplist, file=paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/summary_data/subsample/MIX/",trait,"_UKB_prune_snplist_",i,"_",pop,"_tune_MIX_approx",approx,"_ratio3.00_repeat",rpt,".txt"), 
            row.names=F, col.names=T, quote=F, append=F, sep = "\t")

}

}

}