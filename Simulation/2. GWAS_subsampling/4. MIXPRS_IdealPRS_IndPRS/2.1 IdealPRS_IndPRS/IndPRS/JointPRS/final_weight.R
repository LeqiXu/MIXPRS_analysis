library(data.table)
library(dplyr)

params <- read.table("/gpfs/gibbs/pi/zhao/xz674/simulation/b.method_calculate/IndPRS/JointPRS/params.txt")
colnames(params) <- c("sample1","subpop","sample2","h2","sim_i","rhog","p","ff","chr")
params <- params %>% distinct(across(-chr))

n_trait <- 5
pops <- c("EUR", "EAS", "AFR", "SAS", "AMR")
out_dir <- "/gpfs/gibbs/pi/zhao/xz674/result/sim_result/IndPRS/JointPRS"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

for (i in 1:nrow(params)){
  if (i %% 10 == 0){print(i)}
  param <- params[i, ]
  sample1 <- param$sample1
  subpop <- param$subpop
  sample2 <- param$sample2
  h2 <- param$h2
  sim_i <- param$sim_i
  p <- param$p
  rhog <- param$rhog
  ff <- param$ff
  
  if (subpop == "EUR"){
    pop_pair <- "subEUR_EAS_AFR_SAS_AMR"
  } else if (subpop == "EAS"){
    pop_pair <- "EUR_subEAS_AFR_SAS_AMR"
  } else if (subpop == "AFR"){
    pop_pair <- "EUR_EAS_subAFR_SAS_AMR"
  } else if (subpop == "SAS"){
    pop_pair <- "EUR_EAS_AFR_subSAS_AMR"
  } else if (subpop == "AMR"){
    pop_pair <- "EUR_EAS_AFR_SAS_subAMR"
  }

  out_name <- paste0("sim", sim_i, "_h2", h2, "_p", p, "_rhog", rhog, "_", pop_pair, "_", 
                     sample1, "_", sample2, "_JointPRS_real_fold", ff)
  
  for (chr in 1:22){
    if (chr == 1){
      method_all <- list()
      for (pop in pops){
        method_all[[pop]] <- data.table()
      }
    }
    
    method_trait_chr <- list()
    for (pop in pops){
      method_trait_chr[[pop]] <- fread(paste0(out_dir,"/",out_name,"_",pop,"_pst_eff_a1_b0.5_phiauto_chr",chr,".txt"))
      method_trait_chr[[pop]] <- method_trait_chr[[pop]][, c(2, 4, 6)]
      names(method_trait_chr[[pop]]) = c("rsID", "A1", pop)
      method_all[[pop]] = rbind(method_all[[pop]], method_trait_chr[[pop]])
    }
    
    if (chr == 22){
      for (pop in pops){
        write.table(method_all[[pop]], paste0(out_dir, "/", out_name, "_beta", pop, ".txt"),
                    quote=F,sep='\t',row.names=F,col.names=T)
      }
    }
  }
}
