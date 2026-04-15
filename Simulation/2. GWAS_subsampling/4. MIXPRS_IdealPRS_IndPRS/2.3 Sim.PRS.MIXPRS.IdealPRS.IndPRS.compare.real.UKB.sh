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

    # MIXPRS
    test_MIXPRS = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/MIXPRS/sim",sim_i,"_h2",h2,"_p",p,"_rhog",rhog,"_EUR_EAS_AFR_SAS_AMR_",sample1,"_",sample2,"_JointPRS_SDPRX_prune_snplist_1_non_negative_linear_weights_approxTRUE_MIXPRS_beta_",pop,"_prs_",pop,".sscore"))
    test_MIXPRS = test_MIXPRS[match(test_pheno_id,test_MIXPRS$IID)]

    # IdealPRS
    if (pop == "EUR"){
        test_IdealPRS = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/IdealPRS/sim",sim_i,"_h2",h2,"_p",p,"_rhog",rhog,"_EUR_EAS_AFR_SAS_AMR_",sample1,"_",sample2,"_JointPRS_SDPRX_full_snplist_ukbb_non_negative_linear_weights_approxFALSE_MIXPRS_beta_",pop,"_prs_",pop,".sscore"))
        test_IdealPRS = test_IdealPRS[match(test_pheno_id,test_IdealPRS$IID)]
    } else {
        test_IdealPRS_1kg = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/IdealPRS/sim",sim_i,"_h2",h2,"_p",p,"_rhog",rhog,"_EUR_EAS_AFR_SAS_AMR_",sample1,"_",sample2,"_JointPRS_SDPRX_full_snplist_1kg_non_negative_linear_weights_approxFALSE_MIXPRS_beta_",pop,"_prs_",pop,".sscore"))
        #test_IdealPRS_100K = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/IdealPRS/sim",sim_i,"_h2",h2,"_p",p,"_rhog",rhog,"_EUR_EAS_AFR_SAS_AMR_",sample1,"_",sample2,"_JointPRS_SDPRX_full_snplist_100K_non_negative_linear_weights_approxFALSE_MIXPRS_beta_",pop,"_prs_",pop,".sscore"))
        
        test_IdealPRS = test_IdealPRS_1kg
        test_IdealPRS = test_IdealPRS[match(test_pheno_id,test_IdealPRS$IID)]
    }

    # IndPRS
    test_IndPRS_nnls = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/IndPRS/sim",sim_i,"_h2",h2,"_p",p,"_rhog",rhog,"_EUR_EAS_AFR_SAS_AMR_",sample1,"_",sample2,"_JointPRS_SDPRX_full_snplist_nnls_MIXPRS_beta_",pop,"_prs_",pop,".sscore"))
    test_IndPRS_nnls = test_IndPRS_nnls[match(test_pheno_id,test_IndPRS_nnls$IID)]

    test_IndPRS_ridge = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/IndPRS/sim",sim_i,"_h2",h2,"_p",p,"_rhog",rhog,"_EUR_EAS_AFR_SAS_AMR_",sample1,"_",sample2,"_JointPRS_SDPRX_full_snplist_ridge_MIXPRS_beta_",pop,"_prs_",pop,".sscore"))
    test_IndPRS_ridge = test_IndPRS_ridge[match(test_pheno_id,test_IndPRS_ridge$IID)]

    test_IndPRS_lasso = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/IndPRS/sim",sim_i,"_h2",h2,"_p",p,"_rhog",rhog,"_EUR_EAS_AFR_SAS_AMR_",sample1,"_",sample2,"_JointPRS_SDPRX_full_snplist_lasso_MIXPRS_beta_",pop,"_prs_",pop,".sscore"))
    test_IndPRS_lasso = test_IndPRS_lasso[match(test_pheno_id,test_IndPRS_lasso$IID)]

    test_IndPRS_elasticnet = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/IndPRS/sim",sim_i,"_h2",h2,"_p",p,"_rhog",rhog,"_EUR_EAS_AFR_SAS_AMR_",sample1,"_",sample2,"_JointPRS_SDPRX_full_snplist_elasticnet_MIXPRS_beta_",pop,"_prs_",pop,".sscore"))
    test_IndPRS_elasticnet = test_IndPRS_elasticnet[match(test_pheno_id,test_IndPRS_elasticnet$IID)]

    ## final
    # test result
    test_data = data.table(pheno = scale(test_pheno$pheno),
                           MIXPRS = scale(test_MIXPRS$MIXPRS_AVG),
                           IdealPRS = scale(test_IdealPRS$MIXPRS_AVG),
                           IndPRS_nnls = scale(test_IndPRS_nnls$MIXPRS_AVG),
                           IndPRS_ridge = scale(test_IndPRS_ridge$MIXPRS_AVG),
                           IndPRS_lasso = scale(test_IndPRS_lasso$MIXPRS_AVG),
                           IndPRS_elasticnet = scale(test_IndPRS_elasticnet$MIXPRS_AVG))

    colnames(test_data) = c("pheno","MIXPRS_nnls","IdealPRS_nnls","IndPRS_nnls","IndPRS_ridge","IndPRS_lasso","IndPRS_elasticnet")
    
    sub_prs_table = data.table(n = sim_i, pop = pop, p = p, rhog = rhog, sample1 = sample1, sample2 = sample2, 
                               MIXPRS_nnls = 0, IdealPRS_nnls = 0, IndPRS_nnls = 0, IndPRS_ridge = 0, IndPRS_lasso = 0, IndPRS_elasticnet = 0)

    for (j in 2:7){
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

write.table(prs_table,paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/evaluation/sim_PRS_real_MIXPRS_IdealPRS_IndPRS_r2.csv"),quote=F,sep='\t',row.names=F,col.names=T)
