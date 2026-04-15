## Step2: Obtain final linear combination results
library(data.table)
library(dplyr)

pop1 = "EUR"

linear_weights_table = data.table()

for (pop2 in c("AFR")){
for (trait in c("Height","BMI","SBP","DBP","PLT")){

pop = pop2

scale_pheno = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/pheno_data/",trait,"/",trait,"_scale_",pop,"_doubleid.tsv"))
scale_pheno =scale_pheno[,c(1,3)]
colnames(scale_pheno) = c("eid","pheno")

JointPRS_EUR = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/JointPRS/UKB_",trait,"_JointPRS_real_EUR_EAS_AFR_beta_EUR_prs_",pop2,".sscore"))
JointPRS_EAS = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/JointPRS/UKB_",trait,"_JointPRS_real_EUR_EAS_AFR_beta_EAS_prs_",pop2,".sscore"))
JointPRS_AFR = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/JointPRS/UKB_",trait,"_JointPRS_real_EUR_EAS_AFR_beta_AFR_prs_",pop2,".sscore"))

SDPRX_EUR = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/SDPRX/UKB_",trait,"_SDPRX_real_EUR_",pop2,"_beta_EUR_prs_",pop2,".sscore"))
SDPRX_EAS = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/SDPRX/UKB_",trait,"_SDPRX_real_EUR_EAS_beta_EAS_prs_",pop2,".sscore"))
SDPRX_AFR = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/SDPRX/UKB_",trait,"_SDPRX_real_EUR_AFR_beta_AFR_prs_",pop2,".sscore"))

PRS_table1 = merge(JointPRS_EUR[,c("IID","EUR_AVG")],JointPRS_EAS[,c("IID","EAS_AVG")],by = c("IID"))
PRS_table1 = merge(PRS_table1,JointPRS_AFR[,c("IID","AFR_AVG")],by = c("IID"))
colnames(PRS_table1) = paste0("JointPRS_",colnames(PRS_table1))

PRS_table2 = merge(SDPRX_EUR[,c("IID","EUR_AVG")],SDPRX_EAS[,c("IID","EAS_AVG")],by = c("IID"))
PRS_table2 = merge(PRS_table2,SDPRX_AFR[,c("IID","AFR_AVG")],by = c("IID"))
colnames(PRS_table2) = paste0("SDPRX_",colnames(PRS_table2))

PRS_table = merge(PRS_table1,PRS_table2, by.x = c("JointPRS_IID"), by.y = c("SDPRX_IID"))
pheno_PRS_table = merge(scale_pheno, PRS_table, by.x = c("eid"),by.y = c("JointPRS_IID"))
pheno_PRS_table = pheno_PRS_table[,c(-1)]
pheno_PRS_table = pheno_PRS_table[,c(2:7) := lapply(.SD,scale),.SDcols = c(2:7)]

linear = lm(pheno ~ . + 0, data = pheno_PRS_table)
weights = linear$coefficient
weights[weights<0] = 0
weights = weights / sum(weights)

sub_linear_weights_table = data.table(trait = c(trait), target_pop = c(pop), JointPRS_EUR = weights[1], JointPRS_EAS = weights[2], JointPRS_AFR = weights[3], SDPRX_EUR = weights[4], SDPRX_EAS = weights[5], SDPRX_AFR = weights[6])
linear_weights_table = rbind(linear_weights_table, sub_linear_weights_table)

}
}

write.table(linear_weights_table,paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/Final_weight/diff_cohort/PAGE_simple_linear_weights.txt"),quote=F,sep='\t',row.names=F,col.names=T)
