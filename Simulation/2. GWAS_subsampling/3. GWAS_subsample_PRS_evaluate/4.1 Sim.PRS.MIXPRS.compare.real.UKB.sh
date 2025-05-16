## fivepops
library(data.table)

h2 = 0.4
rhog = 0.8

sample1 = "ukbb"
sample2 = "100K"

prs_table = data.table()

for(pop in c("EUR","EAS","AFR","SAS","AMR")){
for (sim_i in c(1:5)){
for (p in c(0.1, 0.01, 0.001, 5e-04)){

    ## test pheno
    test_pheno = fread(paste0("/gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/pheno_data/",pop,"/test/",pop,"_sim",sim_i,"_h2",h2,"_p",p,"_rhog",rhog,"_10K.phen"),header = FALSE)
    colnames(test_pheno) = c("FID","IID","pheno")
    test_pheno_id = test_pheno$IID

    # MIXPRS_full_FALSE
    test_MIXPRS_full_FALSE = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/MIXPRS/sim",sim_i,"_h2",h2,"_p",p,"_rhog",rhog,"_EUR_EAS_AFR_SAS_AMR_",sample1,"_",sample2,"_JointPRS_SDPRX_full_snplist_linear_weights_approxFALSE_MIXPRS_beta_",pop,"_prs_",pop,".sscore"))
    test_MIXPRS_full_FALSE = test_MIXPRS_full_FALSE[match(test_pheno_id,test_MIXPRS_full_FALSE$IID)]

    # MIXPRS_prune_FALSE
    test_MIXPRS_prune_FALSE = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/MIXPRS/sim",sim_i,"_h2",h2,"_p",p,"_rhog",rhog,"_EUR_EAS_AFR_SAS_AMR_",sample1,"_",sample2,"_JointPRS_SDPRX_prune_snplist_1_linear_weights_approxFALSE_MIXPRS_beta_",pop,"_prs_",pop,".sscore"))
    test_MIXPRS_prune_FALSE = test_MIXPRS_prune_FALSE[match(test_pheno_id,test_MIXPRS_prune_FALSE$IID)]

    # MIXPRS_prune_TRUE
    test_MIXPRS_prune_TRUE = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/MIXPRS/sim",sim_i,"_h2",h2,"_p",p,"_rhog",rhog,"_EUR_EAS_AFR_SAS_AMR_",sample1,"_",sample2,"_JointPRS_SDPRX_prune_snplist_1_linear_weights_approxTRUE_MIXPRS_beta_",pop,"_prs_",pop,".sscore"))
    test_MIXPRS_prune_TRUE = test_MIXPRS_prune_TRUE[match(test_pheno_id,test_MIXPRS_prune_TRUE$IID)]

    ## final
    # test result
    test_data = data.table(pheno = scale(test_pheno$pheno),
                           MIXPRS_full_FALSE = scale(test_MIXPRS_full_FALSE$MIXPRS_AVG),
                           MIXPRS_prune_FALSE = scale(test_MIXPRS_prune_FALSE$MIXPRS_AVG),
                           MIXPRS_prune_TRUE = scale(test_MIXPRS_prune_TRUE$MIXPRS_AVG))
    colnames(test_data) = c("pheno","MIXPRS_full_FALSE","MIXPRS_prune_FALSE","MIXPRS_prune_TRUE")
    
    sub_prs_table = data.table(n = sim_i, pop = pop, p = p, rhog = rhog, sample1 = sample1, sample2 = sample2, 
                               MIXPRS_full_FALSE = 0, MIXPRS_prune_FALSE = 0, MIXPRS_prune_TRUE = 0)

    for (j in 2:4){
        data = data.table(status = test_data$pheno)
        data$prs <- unlist(test_data[,..j])
        linear = lm(status ~ prs, data=data)

        sub_prs_table[1,j+5] = summary(linear)$`r.squared`
    }
    
    prs_table = rbind(prs_table,sub_prs_table)


}
}
}

print(prs_table)

write.table(prs_table,paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/evaluation/sim_PRS_real_MIXPRS_compare_r2.csv"),quote=F,sep='\t',row.names=F,col.names=T)
