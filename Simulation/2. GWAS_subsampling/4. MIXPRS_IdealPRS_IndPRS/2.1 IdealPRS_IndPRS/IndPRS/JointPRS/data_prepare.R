## Step0: Obtain PRScsx format GWAS
## cv IndPRS
library(data.table)

h2 = 0.4
rhog = 0.8

for (pop in c("EUR")){ # ,"EAS","AFR","SAS","AMR"
  sample_size = "ukbb"
  for (sim_i in c(1:5)){
    for (p in c(0.1, 0.01, 0.001, 5e-04)){
      for (ff in c(1:4)){
        
        sumstat_data = fread(paste0("/gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/summary_data/clean/",pop,"_sim",sim_i,"_h2",h2,"_p",p,"_rhog",rhog,"_",sample_size,"_clean_real_fold",ff,".txt"))
        sumstat_data = sumstat_data[,c("SNP","A1","A2","BETA","P")]
        
        write.table(sumstat_data,paste0("/gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/summary_data/PRScsx/",pop,"_sim",sim_i,"_h2",h2,"_p",p,"_rhog",rhog,"_",sample_size,"_PRScsx_real_fold",ff,".txt"),quote=F,sep='\t',row.names=F,col.names=T)
        
      }
    }
  }
}
