# Step1: Data preparation
library(data.table)
library(stringr)
library(dplyr)
library(pROC)

h2 <- 0.4
model_list <- c("linear", "lasso", "ridge", "elasticnet", "nnls")
n_model <- length(model_list)
out_dir <- "/gpfs/gibbs/pi/zhao/xz674/result/sim_result"

for (pop in c("EUR", "EAS", "AFR", "SAS", "AMR")){
  ssize <- ifelse(pop == "EUR", "ukbb", "100K")
  print(paste0("... Train PRS for target pop ", pop, " ..."))
  for (sim in c(1:5)){
    for (p in c(0.001,0.01,0.1,5e-4)){
      for (rhog in c(0.8)){
        # coefs_list <- list()
        # for (model in model_list){
        #   coefs_list[[model]] <- rep(0, 10)
        # }
        for (ff in 1:4){
          setting <- paste0("sim=", sim, ", p=", p, ", rhog=", rhog, ", ff=", ff, ", ssize=", ssize)
          print(paste0("... Setting: ", setting, " ..."))
          
          ## validate phenotype id
          validate_id <- fread(paste0("/gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/info/",
                                      ifelse(pop=="EUR","EUR_",""),"validate_",ssize,"_id_fold",ff,".tsv"))
          validate_id <- validate_id$IID
          pheno_file <- paste0("/gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/pheno_data/",pop,"/discover_validate/",
                               pop,"_sim",sim,"_h2",h2,"_p",p,"_rhog",rhog,"_",ssize,"_all.phen")
          all_pheno <- fread(pheno_file, header=F)
          colnames(all_pheno) <- c("FID", "IID", "pheno")
          pheno_name <- c("IID", "pheno")
          all_pheno <- all_pheno[,..pheno_name]
          validate_pheno <- na.omit(all_pheno[IID %in% validate_id])
          validate_id <- validate_pheno$IID
          
          # PRS list
          PRS_list <- list()
          
          # JointPRS files
          sample1 = "ukbb"
          sample2 = "100K"
          out_dir = "/gpfs/gibbs/pi/zhao/xz674/result/sim_result/IndPRS/JointPRS"
          if (pop == "EUR") {
            out_name <- paste0("sim", sim, "_h2", h2, "_p", p, "_rhog", rhog, "_subEUR_EAS_AFR_SAS_AMR_", sample1, "_", sample2, "_JointPRS_real_fold", ff)
          } else if (pop == "EAS") {
            out_name <- paste0("sim", sim, "_h2", h2, "_p", p, "_rhog", rhog, "_EUR_subEAS_AFR_SAS_AMR_", sample1, "_", sample2, "_JointPRS_real_fold", ff)
          } else if (pop == "AFR") {
            out_name <- paste0("sim", sim, "_h2", h2, "_p", p, "_rhog", rhog, "_EUR_EAS_subAFR_SAS_AMR_", sample1, "_", sample2, "_JointPRS_real_fold", ff)
          } else if (pop == "SAS") {
            out_name <- paste0("sim", sim, "_h2", h2, "_p", p, "_rhog", rhog, "_EUR_EAS_AFR_subSAS_AMR_", sample1, "_", sample2, "_JointPRS_real_fold", ff)
          } else if (pop == "AMR") {
            out_name <- paste0("sim", sim, "_h2", h2, "_p", p, "_rhog", rhog, "_EUR_EAS_AFR_SAS_subAMR_", sample1, "_", sample2, "_JointPRS_real_fold", ff)
          }
          for (pp in c("EUR", "EAS", "AFR", "SAS", "AMR")){
            PRS_list[[paste0("JointPRS_", pp)]] = fread(paste0(out_dir, "/", out_name, "_beta_", pp, "_prs_", pp,".sscore"))[,c("IID","SCORE1_AVG")]
            colnames(PRS_list[[paste0("JointPRS_", pp)]]) = c("IID", paste0("JointPRS_", pp))
          }
          
          # SDPRX files
          out_dir = "/gpfs/gibbs/pi/zhao/xz674/result/sim_result/IndPRS/SDPRX"
          if (pop == "EUR") {
            sample_size <- "ukbb"
            out_name <- paste0("sim", sim, "_h2", h2, "_p", p, "_rhog", rhog, "_subEUR_EAS_", sample1, "_", sample2, "_SDPRX_real_fold", ff)
            PRS_list[["SDPRX_EUR"]] <- fread(paste0(out_dir, "/", out_name, "_beta_", pop, "_prs_", pop,".sscore"))[,c("IID","SCORE1_AVG")]
            colnames(PRS_list[["SDPRX_EUR"]]) = c("IID", "SDPRX_EUR")
            for (pp in c("EAS","AFR","SAS","AMR")){
              out_name <- paste0("sim", sim, "_h2", h2, "_p", p, "_rhog", rhog, "_subEUR_",pp,"_", sample1, "_", sample2, "_SDPRX_real_fold", ff)
              PRS_list[[paste0("SDPRX_", pp)]] <- fread(paste0(out_dir, "/", out_name, "_beta_", pp, "_prs_", pp,".sscore"))[,c("IID","SCORE1_AVG")]
              colnames(PRS_list[[paste0("SDPRX_", pp)]]) = c("IID", paste0("SDPRX_", pp))
            }
          } else {
            sample_size <- "100K"
            out_name <- paste0("sim", sim, "_h2", h2, "_p", p, "_rhog", rhog, "_EUR_sub", pop, "_", sample1, "_", sample2, "_SDPRX_real_fold", ff)
            PRS_list[["SDPRX_EUR"]] <- fread(paste0(out_dir, "/", out_name, "_beta_EUR_prs_EUR.sscore"))[,c("IID","SCORE1_AVG")]
            colnames(PRS_list[["SDPRX_EUR"]]) = c("IID", "SDPRX_EUR")
            for (pp in c("EAS","AFR","SAS","AMR")){
              out_name <- paste0("sim", sim, "_h2", h2, "_p", p, "_rhog", rhog, "_EUR_sub", pp, "_", sample1, "_", sample2, "_SDPRX_real_fold", ff)
              PRS_list[[paste0("SDPRX_", pp)]] <- fread(paste0(out_dir, "/", out_name, "_beta_", pp, "_prs_", pp,".sscore"))[,c("IID","SCORE1_AVG")]
              colnames(PRS_list[[paste0("SDPRX_", pp)]]) = c("IID", paste0("SDPRX_", pp))
            }
          }
          
          methods <- names(PRS_list)
          multi_pop_PRS <- PRS_list[[methods[1]]]
          for (i in 2:length(methods)){
            multi_pop_PRS <- merge(multi_pop_PRS, PRS_list[[methods[i]]], by = c("IID"))
          }
          
          ## match the IID for pheno and PRS fortraining (validate_pheno)
          validate_pheno$pheno = scale(validate_pheno$pheno)
          validate_pheno = validate_pheno[,-1]
          
          multi_pop_PRS_train = multi_pop_PRS[match(validate_id, multi_pop_PRS$IID),-1]
          multi_pop_PRS_train = as.data.table(scale(multi_pop_PRS_train))
          
          # Data for Train
          X_train <- model.matrix( ~ 0 + ., data = multi_pop_PRS_train)
          y_train <- validate_pheno$pheno
          
          # Step2: PRS model training
          out_dir <- "/gpfs/gibbs/pi/zhao/xz674/result/sim_result/IndPRS/Final_weight"
          out_name <- paste0("sim", sim, "_h2", h2, "_p", p, "_rhog", rhog, 
                             "_EUR_EAS_AFR_SAS_AMR_UKB_100K_JointPRS_SDPRX_full_snplist_fold",ff,"_",pop)
          
          for (model in model_list){
            if (model == "lasso"){
              library(glmnet)  # For elastic net
              family_choice <- "gaussian"
              set.seed(123)  # reproducibility
              cvfit <- cv.glmnet(
                x = X_train,
                y = y_train,
                alpha = 1,
                family = family_choice,
                intercept = FALSE  # no intercept in the penalized model
              )
              best_lambda <- cvfit$lambda.min
              coefs <- as.vector(coef(cvfit, s = best_lambda))[-1]
            } else if (model == "ridge"){
              library(glmnet)  # For elastic net
              family_choice <- "gaussian"
              set.seed(123)  # reproducibility
              cvfit <- cv.glmnet(
                x = X_train,
                y = y_train,
                alpha = 0,
                family = family_choice,
                intercept = FALSE  # no intercept in the penalized model
              )
              best_lambda <- cvfit$lambda.min
              coefs <- as.vector(coef(cvfit, s = best_lambda))[-1]
            } else if (model == "elasticnet"){
              library(glmnet)  # For elastic net
              family_choice <- "gaussian"
              set.seed(123)  # reproducibility
              cvfit <- cv.glmnet(
                x = X_train,
                y = y_train,
                alpha = 0.5,
                family = family_choice,
                intercept = FALSE  # no intercept in the penalized model
              )
              best_lambda <- cvfit$lambda.min
              coefs <- as.vector(coef(cvfit, s = best_lambda))[-1]
            } else if (model == "nnls"){
              library(nnls) # For NNLS
              library(glmnet)  # For logistic
              nnfit <- nnls(X_train, y_train)
              coefs <- coef(nnfit)
            } else if (model == "linear"){
              lmfit <- lm(y_train ~ X_train - 1)
              coefs <- as.vector(coef(lmfit))
            } 
            # coefs_list[[model]] <- coefs_list[[model]] + coefs / sum(coefs)
            write.table(t(coefs), file = paste0(out_dir, "/", out_name, "_", model, "_weights.txt"), 
                        row.names = F, col.names = F, sep = "\t")
          }
        }
        # for (model in model_list){
        #   coefs_list[[model]] <- coefs_list[[model]] / ff
        #   write.table(t(coefs_list[[model]]), file = paste0(out_dir, "/", out_name, "_", model, "_weights.txt"), 
        #               row.names = F, col.names = F, sep = "\t")
        # }
      }
    }
  }
}

