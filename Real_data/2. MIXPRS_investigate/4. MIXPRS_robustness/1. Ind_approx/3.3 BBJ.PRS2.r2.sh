# no_val situation
library(data.table)
library(stringr)
library(dplyr)

cov_choice = c("age_recruit","sex",paste0("PC",1:20))
total_covariates = fread("/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/cov_data/adjust.csv")
total_covariates = total_covariates %>% select(all_of(c("eid",cov_choice)))

pop = "EAS"
trait_list = c("WBC","NEU","LYM","MON","EOS","RBC","HCT","MCH","MCV","HB","ALT","ALP","GGT")
n_trait = length(trait_list)
trait_type_list = c(rep("BBJ",13))

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
Trait_pheno$pheno <- scale(Trait_pheno$pheno)
Trait_pheno_id = Trait_pheno$eid

## weights table
weights_table = data.table(GWAS_type = c(GWAS_type), trait = c(trait), target_pop = c(pop), snplist = c(i), approx = c(approx))

## PRS score alignment with trait pheno
Trait_MIXPRS = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/MIXPRS/UKB_",trait,"_MIXPRS_",GWAS_type,"_",i,"_linear_weights_approx",approx,"_prs_",pop,".sscore"))
Trait_MIXPRS = Trait_MIXPRS[match(Trait_pheno_id,Trait_MIXPRS$IID),]

## score alignment with trait pheno
covariates = total_covariates[match(Trait_pheno_id,total_covariates$eid),]
pheno_covariates = cbind(Trait_pheno[,-1],covariates[,-1])
        
# null model in all individuals in UKBB dataset
linear_null = lm(pheno ~ . , data = pheno_covariates)
linear_null_summary = summary(linear_null)
linear_null_res2 = sum(linear_null_summary$residuals^2)
        
# prs comparison            
## weights_table
weights_table$weight_r2 = 0
data = pheno_covariates
data$prs <- scale(Trait_MIXPRS$SCORE1_AVG)
linear = lm(pheno ~ ., data=data)
linear_summary=summary(linear)
linear_summary_res2 = sum(linear_summary$residuals^2)          
weights_table$weight_r2 = 1 - linear_summary_res2/linear_null_res2

prs_table = rbind(prs_table,weights_table)

}
}

write.table(prs_table,paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/evaluation/BBJ_MIX_prune_PRS_r2_ind_approx.csv"),quote=F,sep='\t',row.names=F,col.names=T)