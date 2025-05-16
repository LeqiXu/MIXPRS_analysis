# no_val situation
library(data.table)
library(stringr)
library(dplyr)
library(pROC)

cov_choice = c("age_recruit","sex",paste0("PC",1:20))
total_covariates = fread("/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/cov_data/adjust.csv")
total_covariates = total_covariates %>% select(all_of(c("eid",cov_choice)))

trait_list = c("T2D","BrC")
n_trait = length(trait_list)
trait_type_list = c(rep("Binary_3",2))

prs_table = c()

GWAS_type = "subsample_prune"
i=1

for (pop in c("EAS","AFR")){
for (t in 1:n_trait){
for (approx in c("TRUE","FALSE")){

trait = trait_list[t]
trait_type = trait_type_list[t]

## test phenotype
Trait_pheno <- fread(paste0("/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/pheno_data/",trait,"/",trait,"_",pop,".tsv"))
colnames(Trait_pheno) = c("eid","pheno")
Trait_pheno_id = Trait_pheno$eid

## weights table
weights_table = data.table(GWAS_type = c(GWAS_type), trait = c(trait), target_pop = c(pop), snplist = c(i), approx = c(approx))

## PRS score alignment with trait pheno
Trait_MIXPRS = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/MIXPRS/UKB_",trait,"_MIXPRS_",GWAS_type,"_",i,"_linear_weights_approx",approx,"_prs_",pop,".sscore"))
Trait_MIXPRS = Trait_MIXPRS[match(Trait_pheno_id,Trait_MIXPRS$IID),]

## score alignment with trait pheno
pheno = Trait_pheno[,-1]

# prs comparison            
## weights_table
weights_table$weight_AUC = 0
data = pheno
data$prs <- scale(Trait_MIXPRS$SCORE1_AVG)

glmfit = glm(pheno~prs, data=data,family=binomial(link="logit"))
glmfit_prob = predict(glmfit, type="response")
glmfit_auc = roc(data$pheno, glmfit_prob, quiet=T, plot=F)$auc
weights_table$weight_AUC = glmfit_auc

prs_table = rbind(prs_table,weights_table)

}
}
}

write.table(prs_table,paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/evaluation/Binary_3_MIX_prune_PRS_AUC_ind_approx.csv"),quote=F,sep='\t',row.names=F,col.names=T)


library(data.table)
library(stringr)
library(dplyr)

cov_choice = c("age_recruit","sex",paste0("PC",1:20))
total_covariates = fread("/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/cov_data/adjust.csv")
total_covariates = total_covariates %>% select(all_of(c("eid",cov_choice)))

pop = "EAS"
trait_list = c("CAD","LuC")
n_trait = length(trait_list)
trait_type_list = c(rep("Binary_2",2))

prs_table = c()

GWAS_type = "subsample_prune"
i=1

for (t in 1:n_trait){
for (approx in c("TRUE","FALSE")){

trait = trait_list[t]
trait_type = trait_type_list[t]

## test phenotype
Trait_pheno <- fread(paste0("/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/pheno_data/",trait,"/",trait,"_",pop,".tsv"))
colnames(Trait_pheno) = c("eid","pheno")
Trait_pheno_id = Trait_pheno$eid

## weights table
weights_table = data.table(GWAS_type = c(GWAS_type), trait = c(trait), target_pop = c(pop), snplist = c(i), approx = c(approx))

## PRS score alignment with trait pheno
Trait_MIXPRS = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/MIXPRS/UKB_",trait,"_MIXPRS_",GWAS_type,"_",i,"_linear_weights_approx",approx,"_prs_",pop,".sscore"))
Trait_MIXPRS = Trait_MIXPRS[match(Trait_pheno_id,Trait_MIXPRS$IID),]

## score alignment with trait pheno
pheno = Trait_pheno[,-1]

# prs comparison            
## weights_table
weights_table$weight_AUC = 0
data = pheno
data$prs <- scale(Trait_MIXPRS$SCORE1_AVG)

glmfit = glm(pheno~prs, data=data,family=binomial(link="logit"))
glmfit_prob = predict(glmfit, type="response")
glmfit_auc = roc(data$pheno, glmfit_prob, quiet=T, plot=F)$auc
weights_table$weight_AUC = glmfit_auc

prs_table = rbind(prs_table,weights_table)

}
}

write.table(prs_table,paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/evaluation/Binary_2_MIX_prune_PRS_AUC_ind_approx.csv"),quote=F,sep='\t',row.names=F,col.names=T)