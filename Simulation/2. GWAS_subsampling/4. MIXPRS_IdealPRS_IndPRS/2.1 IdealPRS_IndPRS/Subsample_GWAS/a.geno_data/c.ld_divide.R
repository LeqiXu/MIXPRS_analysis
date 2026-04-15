## 1. Divide LD blocks for each population
library(data.table)

print("Start!")

for (pop in c("EUR", "AFR", "EAS", "SAS", "AMR")){
  if (pop == "EUR"){
    blk = read.table("/gpfs/gibbs/pi/zhao/xz674/data/ref_data/EUR_block.txt", header=T)
  }else if (pop %in% c("SAS", "EAS")){
    blk = read.table("/gpfs/gibbs/pi/zhao/xz674/data/ref_data/EAS_block.txt", header=T)
  }else if(pop %in% c("AFR", "AMR")){
    blk = read.table("/gpfs/gibbs/pi/zhao/xz674/data/ref_data/AFR_block.txt", header=T)
  }
  
  # for (geno in c("100K", "1kg")){
  for (geno in c("100K")){
    print(paste0("Start: pop=", pop, ", geno=", geno))
    
    # # discover_validate
    # print("discover_validate")
    # data = read.table(paste0("/gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/",pop,"/discover_validate/",pop,"_",geno,"_all.bim"))
    # colnames(data) = c("chr","snpid","nonthing","position","minor","major")
    
    # size = rep(0, nrow(blk))
    # chrom = rep(0, nrow(blk))
    # for (i in 1:nrow(blk)) {
    #   chr = as.numeric(sub("chr","",blk[i,1]))
    #   idx = which(data$position >= as.numeric(blk[i,2]) & data$position < as.numeric(blk[i,3]) & data$chr == chr)
    #   if (length(idx) > 0) {
    #     write.table(data[idx,2], file=paste0("/gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/ref_data/ldblk_snplist/discover_validate/",geno,"/",pop,"/",i), 
    #                 row.names=F, col.names=F, quote=F, append=F)
    #   }
    #   size[i] = length(idx)
    #   chrom[i] = chr
    # }
    
    # write.table(size, file=paste0("/gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/ref_data/ldblk_info/discover_validate/",pop,"_",geno,"_blk_size.txt"), 
    #             row.names=F, col.names=F, quote=F, append=F)
    # write.table(chrom, file=paste0("/gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/ref_data/ldblk_info/discover_validate/",pop,"_",geno,"_blk_chr.txt"), 
    #             row.names=F, col.names=F, quote=F, append=F)
    
    # discover
    for (ff in 1:4){
      print(paste0("discover_fold",ff))
      data = read.table(paste0("/gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/",pop,"/discover/",pop,"_",geno,"_fold",ff,".bim"))
      colnames(data) = c("chr","snpid","nonthing","position","minor","major")
      
      size = rep(0, nrow(blk))
      chrom = rep(0, nrow(blk))
      for (i in 1:nrow(blk)) {
        chr = as.numeric(sub("chr","",blk[i,1]))
        idx = which(data$position >= as.numeric(blk[i,2]) & data$position < as.numeric(blk[i,3]) & data$chr == chr)
        if (length(idx) > 0) {
          write.table(data[idx,2], file=paste0("/gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/ref_data/ldblk_snplist/discover/",geno,"_fold",ff,"/",pop,"/",i), 
                      row.names=F, col.names=F, quote=F, append=F)
        }
        size[i] = length(idx)
        chrom[i] = chr
      }
      
      write.table(size, file=paste0("/gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/ref_data/ldblk_info/discover/",pop,"_",geno,"_blk_size_fold",ff,".txt"), 
                  row.names=F, col.names=F, quote=F, append=F)
      write.table(chrom, file=paste0("/gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/ref_data/ldblk_info/discover/",pop,"_",geno,"_blk_chr_fold",ff,".txt"), 
                  row.names=F, col.names=F, quote=F, append=F)
    }
  }
}

print("Finish!")