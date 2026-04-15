library(data.table)
library(dplyr)

params <- read.table("/gpfs/gibbs/pi/zhao/xz674/simulation/b.method_calculate/IndPRS/SDPRX/params.txt")
colnames(params) <- c("pop1","sample1","pop2","sample2","subpop","h2","sim_i","rhog","p","ff","chr")
params <- params %>% distinct(across(-chr))

n_pop <- 2
out_dir <- "/gpfs/gibbs/pi/zhao/xz674/result/sim_result/IndPRS/SDPRX"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

for (i in 1:nrow(params)){
  if (i %% 10 == 0){print(i)}
  param <- params[i, ]
  pop1 <- param$pop1
  sample1 <- param$sample1
  pop2 <- param$pop2
  sample2 <- param$sample2
  subpop <- param$subpop
  h2 <- param$h2
  sim_i <- param$sim_i
  p <- param$p
  rhog <- param$rhog
  ff <- param$ff
  
  if (subpop == "EUR"){
    out_name <- paste0("sim",sim_i,"_h2",h2,"_p",p,"_rhog",rhog,"_sub",pop1,"_",pop2,"_",sample1,"_",sample2,
                       "_SDPRX_real_fold",ff)
  } else{
    out_name <- paste0("sim",sim_i,"_h2",h2,"_p",p,"_rhog",rhog,"_",pop1,"_sub",pop2,"_",sample1,"_",sample2,
                       "_SDPRX_real_fold",ff)
  }
  
  for (chr in 1:22){
    if (chr == 1){
      method_all <- list()
      for (pp in 1:n_pop){
        method_all[[pp]] <- data.table()
      }
    }
    
    method_trait_chr <- list()
    for (pp in 1:n_pop){
      method_trait_chr[[pp]] <- fread(paste0(out_dir,"/",out_name,"_chr",chr,"_",pp,".txt"))
      pop = ifelse(pp == 1, pop1, pop2)
      names(method_trait_chr[[pp]]) = c("rsID", "A1", pop)
      method_all[[pp]] = rbind(method_all[[pp]], method_trait_chr[[pp]])
    }
    
    if (chr == 22){
      for (pp in 1:n_pop){
        pop = ifelse(pp == 1, pop1, pop2)
        write.table(method_all[[pp]], paste0(out_dir, "/", out_name, "_beta", pop, ".txt"),
                    quote=F,sep='\t',row.names=F,col.names=T)
      }
    }
  }
}

