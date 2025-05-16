### Step1: Organize the result
library(data.table)

h2 = 0.4
rhog = 0.8
i=1

corr_table = data.table()

for (pop in c("EUR","EAS","AFR","SAS","AMR")){

if (pop == "EUR"){
sample_size = "ukbb"
}

if (pop != "EUR"){
sample_size = "100K"
}

for (sim_i in c(1:5)){
for (repeat_i in c(1:4)){
for (p in c(0.1, 0.01, 0.001, 5e-04)){

print(paste0(pop, "_sim",sim_i,"_repeat",repeat_i,"_p",p))

combined_residual_full_snplist_train = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/residual/",pop,"_sim",sim_i,"_h2",h2,"_p",p,"_rhog",rhog,"_",sample_size,"_full_snplist_",pop,"_train_GWAS_residual_approxFALSE_ratio3.00_repeat",repeat_i,".txt"))
combined_residual_full_snplist_tune = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/residual/",pop,"_sim",sim_i,"_h2",h2,"_p",p,"_rhog",rhog,"_",sample_size,"_full_snplist_",pop,"_tune_GWAS_residual_approxFALSE_ratio3.00_repeat",repeat_i,".txt"))
combined_residual_full_snplist = merge(combined_residual_full_snplist_train[,c("SNP","Residual")],combined_residual_full_snplist_tune[,c("SNP","Residual")], by = c("SNP"))
colnames(combined_residual_full_snplist) = c("SNP","train_residual","tune_residual")
combined_residual_full_snplist <- combined_residual_full_snplist[is.finite(combined_residual_full_snplist$train_residual) & is.finite(combined_residual_full_snplist$tune_residual),]

print("full_snplist")
print(var(combined_residual_full_snplist$train_residual))
print(var(combined_residual_full_snplist$tune_residual))

sub_corr_table = data.table(pop = c(pop), sample_size = c(sample_size), sim_i = c(sim_i), repeat_i = c(repeat_i), p = c(p), type = c("full_snplist"), corr = c(cor(combined_residual_full_snplist$train_residual,combined_residual_full_snplist$tune_residual)))
corr_table = rbind(corr_table,sub_corr_table)

combined_residual_prune_snplist_train = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/residual/",pop,"_sim",sim_i,"_h2",h2,"_p",p,"_rhog",rhog,"_",sample_size,"_prune_snplist_",i,"_",pop,"_train_GWAS_residual_approxFALSE_ratio3.00_repeat",repeat_i,".txt"))
combined_residual_prune_snplist_tune = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/residual/",pop,"_sim",sim_i,"_h2",h2,"_p",p,"_rhog",rhog,"_",sample_size,"_prune_snplist_",i,"_",pop,"_tune_GWAS_residual_approxFALSE_ratio3.00_repeat",repeat_i,".txt"))
combined_residual_prune_snplist = merge(combined_residual_prune_snplist_train[,c("SNP","Residual")],combined_residual_prune_snplist_tune[,c("SNP","Residual")], by = c("SNP"))
colnames(combined_residual_prune_snplist) = c("SNP","train_residual","tune_residual")
combined_residual_prune_snplist <- combined_residual_prune_snplist[is.finite(combined_residual_prune_snplist$train_residual) & is.finite(combined_residual_prune_snplist$tune_residual),]

print("prune_snplist")
print(var(combined_residual_prune_snplist$train_residual))
print(var(combined_residual_prune_snplist$tune_residual))

sub_corr_table = data.table(pop = c(pop), sample_size = c(sample_size), sim_i = c(sim_i), repeat_i = c(repeat_i), p = c(p), type = c("prune_snplist"), corr = c(cor(combined_residual_prune_snplist$train_residual,combined_residual_prune_snplist$tune_residual)))
corr_table = rbind(corr_table,sub_corr_table)

combined_residual_prune_snplist_ind_approx_train = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/residual/",pop,"_sim",sim_i,"_h2",h2,"_p",p,"_rhog",rhog,"_",sample_size,"_prune_snplist_",i,"_",pop,"_train_GWAS_residual_approxTRUE_ratio3.00_repeat",repeat_i,".txt"))
combined_residual_prune_snplist_ind_approx_tune = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/residual/",pop,"_sim",sim_i,"_h2",h2,"_p",p,"_rhog",rhog,"_",sample_size,"_prune_snplist_",i,"_",pop,"_tune_GWAS_residual_approxTRUE_ratio3.00_repeat",repeat_i,".txt"))
combined_residual_prune_snplist_ind_approx = merge(combined_residual_prune_snplist_ind_approx_train[,c("SNP","Residual")],combined_residual_prune_snplist_ind_approx_tune[,c("SNP","Residual")], by = c("SNP"))
colnames(combined_residual_prune_snplist_ind_approx) = c("SNP","train_residual","tune_residual")
combined_residual_prune_snplist_ind_approx <- combined_residual_prune_snplist_ind_approx[is.finite(combined_residual_prune_snplist_ind_approx$train_residual) & is.finite(combined_residual_prune_snplist_ind_approx$tune_residual),]

print("prune_snplist_ind_approx")
print(var(combined_residual_prune_snplist_ind_approx$train_residual))
print(var(combined_residual_prune_snplist_ind_approx$tune_residual))

sub_corr_table = data.table(pop = c(pop), sample_size = c(sample_size), sim_i = c(sim_i), repeat_i = c(repeat_i), p = c(p), type = c("prune_snplist_ind_approx"), corr = c(cor(combined_residual_prune_snplist_ind_approx$train_residual,combined_residual_prune_snplist_ind_approx$tune_residual)))
corr_table = rbind(corr_table,sub_corr_table)

}
}
}
}

write.table(corr_table,paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/residual/simulation_residual_corr.txt"),quote=F,sep='\t',row.names=F,col.names=T)
