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

    # JointPRS-auto
    test_JointPRS_auto = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/JointPRS/sim",sim_i,"_h2",h2,"_p",p,"_rhog",rhog,"_",sample1,"_",sample2,"_JointPRS_real_EUR_EAS_AFR_SAS_AMR_beta_",pop,"_prs_",pop,".sscore"))
    test_JointPRS_auto = test_JointPRS_auto[match(test_pheno_id,test_JointPRS_auto$IID),]

    # SDPRX
    if (pop == "EUR"){
        test_SDPRX = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/SDPRX/sim",sim_i,"_h2",h2,"_p",p,"_rhog",rhog,"_",sample1,"_",sample2,"_SDPRX_real_EUR_EAS_beta_EUR_prs_EUR.sscore"))
    } else {
        test_SDPRX = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/SDPRX/sim",sim_i,"_h2",h2,"_p",p,"_rhog",rhog,"_",sample1,"_",sample2,"_SDPRX_real_EUR_",pop,"_beta_",pop,"_prs_",pop,".sscore"))
    }
    test_SDPRX = test_SDPRX[match(test_pheno_id,test_SDPRX$IID),]

    # MIXPRS_prune_TRUE
    test_MIXPRS_prune_TRUE = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/MIXPRS/sim",sim_i,"_h2",h2,"_p",p,"_rhog",rhog,"_EUR_EAS_AFR_SAS_AMR_",sample1,"_",sample2,"_JointPRS_SDPRX_prune_snplist_1_non_negative_linear_weights_approxTRUE_MIXPRS_beta_",pop,"_prs_",pop,".sscore"))
    test_MIXPRS_prune_TRUE = test_MIXPRS_prune_TRUE[match(test_pheno_id,test_MIXPRS_prune_TRUE$IID)]

    ## final
    # test result
    test_data = data.table(pheno = scale(test_pheno$pheno),
                           JointPRS_auto = unlist(scale(test_JointPRS_auto[,5])),
                           SDPRX = unlist(scale(test_SDPRX[,5])),
                           MIXPRS_prune_TRUE = scale(test_MIXPRS_prune_TRUE$MIXPRS_AVG))
    colnames(test_data) = c("pheno","JointPRS_auto","SDPRX","MIXPRS_prune_TRUE")
    
    sub_prs_table = data.table(n = sim_i, pop = pop, p = p, rhog = rhog, sample1 = sample1, sample2 = sample2, 
                               JointPRS_auto = 0, SDPRX = 0, MIXPRS_prune_TRUE = 0)

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

write.table(prs_table,paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/evaluation/sim_PRS_real_MIXPRS_benefit_r2.csv"),quote=F,sep='\t',row.names=F,col.names=T)
