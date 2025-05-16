# Step1: Organize beta by chr pop for each param in each trait
library(data.table)

for (trait in c("HDL","LDL","TC","logTG")){
for (pop in c("EUR","EAS","AFR","SAS")){

PUMAS_all <- fread(paste0("/vast/palmer/pi/zhao/yd357/pumas_evaluation_output_filtered/",trait,"/",pop,"/",trait,"_",pop,"_inter.ensemble.weights.txt"), sep="\t", header=TRUE, fill=TRUE)
PUMAS_all <- PUMAS_all[,c("SNP","A1","EN")]

write.table(PUMAS_all,paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/PUMAS_EN/",trait,"_PUMAS_EN_EUR_EAS_AFR_SAS_beta_",pop,".txt"),quote=F,sep='\t',row.names=F,col.names=T)

}
}

