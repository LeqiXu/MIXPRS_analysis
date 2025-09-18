## 1. Auto PRS methods test and compare
## fivepops
library(data.table)

h2 = 0.4
rho_list = c(0.8) 
p_list = c(0.001,0.01,5e-04,0.1)

sample1 = "UKB"
sample2_list = c("25K","90K")

pop_list = c("EUR","EAS","AFR","SAS","AMR")
prs_table = data.table()

for(pop in c("EAS","AFR","SAS","AMR")){
for(sim_i in c(1:5)){
for (p in p_list){
for (rhog in rho_list){
for (sample2 in sample2_list){

    if (sample2 == "25K"){
        sample2_train = "15K"
        sample2_val = "10K"
    }

    if (sample2 == "90K"){
        sample2_train = "80K"
        sample2_val = "10K"
    }

    sub_prs_table = data.table(n = c(sim_i), pop = c(pop), p = c(p), rhog = c(rhog), sample2 = c(sample2),
                               MIXPRS_auto_5 = c(0), SDPRX_auto_2 =c(0), XPASS_auto_2 =c(0), JointPRS_tune_5 =c(0), PRScsx_tune_5 =c(0), PROSPER_tune_5 =c(0), MUSSEL_tune_5 =c(0), BridgePRS_tune_2 =c(0))

    ## test pheno
    test_pheno = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/data/sim_data/pheno_data/",pop,"/test/",pop,"_sim",sim_i,"_p",p,"_rho",rhog,"_10K_doubleidname.tsv"))
    test_pheno_id = test_pheno$IID

    ## auto methods
    # MIXPRS_5
    type = "prune_snplist_1"; approx = "TRUE"
    test_MIXPRS_5 = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/MIXPRS/test_sim",sim_i,"_h2",h2,"_p",p,"_rhog",rhog,"_EUR_EAS_AFR_SAS_AMR_",sample1,"_",sample2,"_JointPRS_SDPRX_",type,"_non_negative_linear_weights_approx",approx,"_",pop,"_MIXPRS_prs_",pop,".sscore"))
    test_MIXPRS_5 = test_MIXPRS_5[match(test_pheno_id,test_MIXPRS_5$IID),]

    # SDPRX_2
    test_SDPRX_2 = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/SDPRX/test_sim",sim_i,"_h2",h2,"_p",p,"_rhog",rhog,"_",sample1,"_",sample2,"_SDPRX_real_EUR_",pop,"_beta_",pop,"_prs_",pop,".sscore"))
    test_SDPRX_2 = test_SDPRX_2[match(test_pheno_id,test_SDPRX_2$IID),]

    # XPASS_2
    test_XPASS_2 = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/result/sim_result/XPASS/sim",sim_i,"_p",p,"_rho",rhog,"/test_sim",sim_i,"_p",p,"_rho",rhog,"_",sample1,"_",sample2,"_XPASS_real_EUR_",pop,"_prs_",pop,".sscore"))
    test_XPASS_2 = test_XPASS_2[match(test_pheno_id,test_XPASS_2$IID),]

    ## tune methods
    R2_sub = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/result/sim_result/Final_weight/auto/JointPRS/sim",sim_i,"_p",p,"_rho",rhog,"_UKB_",sample2_train,"_",sample2_val,"_JointPRS_real_EUR_EAS_AFR_SAS_AMR_r2_",pop,".txt"))
    R2_full = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/result/sim_result/Final_weight/tuning/JointPRS_tune/sim",sim_i,"_p",p,"_rho",rhog,"_UKB_",sample2_train,"_",sample2_val,"_JointPRS_real_EUR_EAS_AFR_SAS_AMR_r2_",pop,".txt"))
    p_value = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/result/sim_result/Final_weight/auto/JointPRS/sim",sim_i,"_p",p,"_rho",rhog,"_UKB_",sample2_train,"_",sample2_val,"_JointPRS_real_EUR_EAS_AFR_SAS_AMR_pvalue_",pop,".txt"))

    if (R2_sub > 0.01 && (R2_full - R2_sub > 0) && p_value < 0.05){
    JointPRS_weight = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/result/sim_result/Final_weight/tuning/JointPRS_tune/sim",sim_i,"_p",p,"_rho",rhog,"_",sample1,"_",sample2_train,"_",sample2_val,"_JointPRS_real_EUR_EAS_AFR_SAS_AMR_weight_",pop,".txt"))

    test_JointPRS_5_EUR = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/JointPRS/test_sim",sim_i,"_h2",h2,"_p",p,"_rhog",rhog,"_",sample1,"_",sample2,"_JointPRS_real_EUR_EAS_AFR_SAS_AMR_beta_EUR_prs_",pop,".sscore"))
    test_JointPRS_5_EUR = test_JointPRS_5_EUR[match(test_pheno_id,test_JointPRS_5_EUR$IID),]

    test_JointPRS_5_EAS = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/JointPRS/test_sim",sim_i,"_h2",h2,"_p",p,"_rhog",rhog,"_",sample1,"_",sample2,"_JointPRS_real_EUR_EAS_AFR_SAS_AMR_beta_EAS_prs_",pop,".sscore"))
    test_JointPRS_5_EAS = test_JointPRS_5_EAS[match(test_pheno_id,test_JointPRS_5_EAS$IID),]
                                                    
    test_JointPRS_5_AFR = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/JointPRS/test_sim",sim_i,"_h2",h2,"_p",p,"_rhog",rhog,"_",sample1,"_",sample2,"_JointPRS_real_EUR_EAS_AFR_SAS_AMR_beta_AFR_prs_",pop,".sscore"))
    test_JointPRS_5_AFR = test_JointPRS_5_AFR[match(test_pheno_id,test_JointPRS_5_AFR$IID),]

    test_JointPRS_5_SAS = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/JointPRS/test_sim",sim_i,"_h2",h2,"_p",p,"_rhog",rhog,"_",sample1,"_",sample2,"_JointPRS_real_EUR_EAS_AFR_SAS_AMR_beta_SAS_prs_",pop,".sscore"))
    test_JointPRS_5_SAS = test_JointPRS_5_SAS[match(test_pheno_id,test_JointPRS_5_SAS$IID),]
                                                    
    test_JointPRS_5_AMR = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/JointPRS/test_sim",sim_i,"_h2",h2,"_p",p,"_rhog",rhog,"_",sample1,"_",sample2,"_JointPRS_real_EUR_EAS_AFR_SAS_AMR_beta_AMR_prs_",pop,".sscore"))
    test_JointPRS_5_AMR = test_JointPRS_5_AMR[match(test_pheno_id,test_JointPRS_5_AMR$IID),]

    test_JointPRS_5 = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/JointPRS/test_sim",sim_i,"_h2",h2,"_p",p,"_rhog",rhog,"_",sample1,"_",sample2,"_JointPRS_real_EUR_EAS_AFR_SAS_AMR_beta_",pop,"_prs_",pop,".sscore"))
    test_JointPRS_5 = test_JointPRS_5[,c("IID")]
    test_JointPRS_5$SCORE1_AVG = JointPRS_weight$EUR * scale(test_JointPRS_5_EUR$SCORE1_AVG) + JointPRS_weight$EAS * scale(test_JointPRS_5_EAS$SCORE1_AVG) + JointPRS_weight$AFR * scale(test_JointPRS_5_AFR$SCORE1_AVG) + JointPRS_weight$SAS * scale(test_JointPRS_5_SAS$SCORE1_AVG) + JointPRS_weight$AMR * scale(test_JointPRS_5_AMR$SCORE1_AVG)
    } else {
    test_JointPRS_5 = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/JointPRS/test_sim",sim_i,"_h2",h2,"_p",p,"_rhog",rhog,"_",sample1,"_",sample2,"_JointPRS_real_EUR_EAS_AFR_SAS_AMR_beta_",pop,"_prs_",pop,".sscore"))
    test_JointPRS_5 = test_JointPRS_5[match(test_pheno_id,test_JointPRS_5$IID),]
    }

    # PRScsx_5
    test_PRScsx_5 = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/result/sim_result/PRScsx/sim",sim_i,"_p",p,"_rho",rhog,"/test_sim",sim_i,"_p",p,"_rho",rhog,"_",sample1,"_",sample2_train,"_",sample2_val,"_PRScsx_real_EUR_EAS_AFR_SAS_AMR_prs_",pop,".sscore"))
    test_PRScsx_5 = test_PRScsx_5[match(test_pheno_id,test_PRScsx_5$IID),]
    PRScsx_weight = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/result/sim_result/Final_weight/tuning/PRScsx/sim",sim_i,"_p",p,"_rho",rhog,"_",sample1,"_",sample2_train,"_",sample2_val,"_PRScsx_real_EUR_EAS_AFR_SAS_AMR_weight_",pop,".txt"))
    test_PRScsx_5$SCORE1_AVG = unlist(PRScsx_weight[1,1]) * unlist(scale(test_PRScsx_5[,5])) + unlist(PRScsx_weight[1,2]) * unlist(scale(test_PRScsx_5[,6])) + unlist(PRScsx_weight[1,3]) * unlist(scale(test_PRScsx_5[,7])) + unlist(PRScsx_weight[1,4]) * unlist(scale(test_PRScsx_5[,8])) + unlist(PRScsx_weight[1,5]) * unlist(scale(test_PRScsx_5[,9]))

    # PROSPER_5
    test_PROSPER_5 = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/result/sim_result/PROSPER/sim",sim_i,"_p",p,"_rho",rhog,"/test_sim",sim_i,"_p",p,"_rho",rhog,"_",sample1,"_",sample2_train,"_",sample2_val,"_PROSPER_update_real_EUR_EAS_AFR_SAS_AMR_prs_",pop,".sscore"))
    test_PROSPER_5 = test_PROSPER_5[match(test_pheno_id,test_PROSPER_5$IID),]

    # MUSSEL_5
    test_MUSSEL_5 = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/result/sim_result/MUSSEL/sim",sim_i,"_p",p,"_rho",rhog,"/test_sim",sim_i,"_p",p,"_rho",rhog,"_",sample1,"_",sample2_train,"_",sample2_val,"_MUSSEL_real_EUR_EAS_AFR_SAS_AMR_prs_",pop,".sscore"))
    test_MUSSEL_5 = test_MUSSEL_5[match(test_pheno_id,test_MUSSEL_5$IID),]

    # BridgePRS_2
    test_BridgePRS_2 = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/result/sim_result/BridgePRS/sim",sim_i,"_p",p,"_rho",rhog,"/test_sim",sim_i,"_p",p,"_rho",rhog,"_",sample1,"_",sample2_train,"_",sample2_val,"_BridgePRS_real_EUR_",pop,"_prs_",pop,".sscore"))
    test_BridgePRS_2 = test_BridgePRS_2[match(test_pheno_id,test_BridgePRS_2$IID),]

    ## final
    # test result
    test_data = data.table(pheno = scale(test_pheno$pheno),
                           MIXPRS_auto_5 = scale(test_MIXPRS_5$SCORE1_AVG),
                           SDPRX_auto_2 = scale(test_SDPRX_2$SCORE1_AVG),
                           XPASS_auto_2 = scale(test_XPASS_2$SCORE1_AVG),
                           JointPRS_tune_5 = scale(test_JointPRS_5$SCORE1_AVG), 
                           PRScsx_tune_5 = scale(test_PRScsx_5$SCORE1_AVG), 
                           PROSPER_tune_5 = scale(test_PROSPER_5$SCORE1_AVG), 
                           MUSSEL_tune_5 = scale(test_MUSSEL_5$SCORE1_AVG),
                           BridgePRS_tune_2 = scale(test_BridgePRS_2$SCORE1_AVG))

    colnames(test_data) = c("pheno","MIXPRS_auto_5","SDPRX_auto_2","XPASS_auto_2","JointPRS_tune_5","PRScsx_tune_5","PROSPER_tune_5","MUSSEL_tune_5","BridgePRS_tune_2")
    for (j in 2:9){
        data = data.table(status = test_data$pheno)
        data$prs <- unlist(test_data[,..j])
        linear = lm(status ~ prs, data=data)
        sub_prs_table[1,j+4] = summary(linear)$`r.squared`
    }
    prs_table = rbind(sub_prs_table,prs_table)

}
}
}
}
}

write.table(prs_table,paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/evaluation/sim_PRS_real_tune_benchmark_r2.csv"),quote=F,sep='\t',row.names=F,col.names=T)