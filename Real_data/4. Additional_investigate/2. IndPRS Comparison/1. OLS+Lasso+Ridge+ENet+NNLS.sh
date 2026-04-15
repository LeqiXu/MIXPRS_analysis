library(data.table)
library(dplyr)
library(caret)
library(glmnet)
library(nnls)

configs <- list(
  GLGC = list(
    traits = c("HDL", "LDL", "TC", "logTG"),
    pops = c("EAS", "AFR", "SAS", "AMR"),
    pop_suffix = "EUR_EAS_AFR_SAS_AMR", 
    weight_pops = c("EUR", "EAS", "AFR", "SAS", "AMR")
  ),
  PAGE = list(
    traits = c("Height", "BMI", "SBP", "DBP", "PLT"),
    pops = c("EAS", "AFR"),
    pop_suffix = "EUR_EAS_AFR",
    weight_pops = c("EUR", "EAS", "AFR")
  ),
  BBJ = list(
    traits = c("WBC", "NEU", "LYM", "MON", "EOS", "RBC", "HCT", "MCH", "MCV", "HB", "ALT", "ALP", "GGT"),
    pops = c("EAS"),
    pop_suffix = "EUR_EAS",
    weight_pops = c("EUR", "EAS")
  )
)

cov_choice <- c("age_recruit", "sex", paste0("PC", 1:20))
total_covariates <- fread("/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/cov_data/adjust.csv")
total_covariates <- total_covariates %>% select(all_of(c("eid", cov_choice)))

results_df <- data.frame()

for (group in names(configs)) {
  
  config <- configs[[group]]
  traits <- config$traits
  pops   <- config$pops
  suffix <- config$pop_suffix
  weight_pops <- config$weight_pops
  
  for (trait in traits) {
    for (pop in pops) {
      
      cat(sprintf("\nProcessing Group: %s | Trait: %s | Target Pop: %s\n", group, trait, pop))
      
      pheno_file <- paste0("/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/pheno_data/", trait, "/", trait, "_", pop, ".tsv")
      pheno_df <- fread(pheno_file)
      colnames(pheno_df)[2] <- "phenotype" 
      
      ## score alignment with trait pheno (Merge Pheno + Covariates)
      model_data <- inner_join(pheno_df, total_covariates, by = "eid")
      
      # MIXPRS
      mix_df <- fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/MIXPRS/UKB_", trait, "_MIXPRS_subsample_prune_prs_", pop, ".sscore"))
      mix_df <- mix_df %>% select(IID, SCORE1_AVG) %>% rename(eid = IID, MIXPRS = SCORE1_AVG)
      model_data <- inner_join(model_data, mix_df, by = "eid")
      
      # Add scores components to IndPRS (JointPRS + SDPRX for all weight pops)
      indprs_col_names <- c()
      
      # Loop through WEIGHT POP (EUR, EAS, etc.)
      for (w_pop in weight_pops) {
        
        # pattern for JointPRS
        jp_df <- fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/JointPRS/UKB_", 
                              trait, "_JointPRS_real_", suffix, "_beta_", w_pop, "_prs_", pop, ".sscore"))
        target_col <- paste0(w_pop, "_AVG")
        new_col_name <- paste0("JointPRS_", w_pop)
        jp_df <- jp_df %>% select(IID, all_of(target_col)) %>% rename(eid = IID, !!new_col_name := !!target_col)
        model_data <- inner_join(model_data, jp_df, by = "eid")
        indprs_col_names <- c(indprs_col_names, new_col_name)
        
        # pattern for SDPRX        
        if (w_pop == "EUR") {
          sd_file_name <- paste0("UKB_", trait, "_SDPRX_real_EUR_", pop, "_beta_EUR_prs_", pop, ".sscore")
        } else {
          sd_file_name <- paste0("UKB_", trait, "_SDPRX_real_EUR_", w_pop, "_beta_", w_pop, "_prs_", pop, ".sscore")
        }
        sd_df <- fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/SDPRX/", sd_file_name))
        target_col <- paste0(w_pop, "_AVG")
        new_col_name <- paste0("SDPRX_", w_pop)
        sd_df <- sd_df %>% select(IID, all_of(target_col)) %>% rename(eid = IID, !!new_col_name := !!target_col)
        model_data <- inner_join(model_data, sd_df, by = "eid")
        indprs_col_names <- c(indprs_col_names, new_col_name)
      }
      

      n_samples <- nrow(model_data)
     

      # Prepare Data for GLMNET (Matrix Format)
      formula_str <- paste("~", paste(c(cov_choice, indprs_col_names), collapse = " + "))
      
      X_full <- model.matrix(as.formula(formula_str), data = model_data)[, -1]
      y_full <- model_data$phenotype
      
      var_names <- colnames(X_full)
      p_factor <- ifelse(grepl("JointPRS_|SDPRX_", var_names), 1, 0) # p_factor: 1 for PRS (penalize), 0 for Covariates
      
      set.seed(1234)
      folds <- createFolds(y_full, k = 4, list = TRUE, returnTrain = FALSE) # 4 Fold CV
      mix_r2_folds   <- numeric(4) # MIXPRS
      ind_r2_folds   <- numeric(4) # OLS
      lasso_r2_folds <- numeric(4) # Lasso
      ridge_r2_folds <- numeric(4) # Ridge
      enet_r2_folds  <- numeric(4) # Elastic Net
      nnls_r2_folds  <- numeric(4) # NNLS
      
for (i in 1:4) {
  test_idx   <- folds[[i]]
  train_idx  <- setdiff(1:n_samples, test_idx)
  
  # Split Dataframes (for lm and nnls), standardizations
        train_df <- model_data[train_idx, ]
        test_df  <- model_data[test_idx, ]
        prs_means <- colMeans(train_df[, ..indprs_col_names])
        prs_sds <- apply(train_df[, ..indprs_col_names], 2, sd)
        for (col in indprs_col_names) { # indprs_col_names is the PRS models we gathered for combining
          train_df[[col]] <- (train_df[[col]] - prs_means[[col]]) / prs_sds[[col]] # Standardization
          test_df[[col]] <- (test_df[[col]] - prs_means[[col]]) / prs_sds[[col]] # Standardization
        }   
        mixprs_mean <- mean(train_df$MIXPRS)
        mixprs_sd <- sd(train_df$MIXPRS)
        train_df$MIXPRS <- (train_df$MIXPRS - mixprs_mean) / mixprs_sd  # Standardization
        test_df$MIXPRS <- (test_df$MIXPRS - mixprs_mean) / mixprs_sd  # Standardization

        x_train <- model.matrix(as.formula(formula_str), data = train_df)[, -1]
        y_train <- train_df$phenotype
        x_test  <- model.matrix(as.formula(formula_str), data = test_df)[, -1]
        y_test  <- test_df$phenotype
        
  
  # 1. Null 
  null_mod  <- lm(as.formula(paste("phenotype ~", paste(cov_choice, collapse="+"))), data = train_df)
  pred_null <- predict(null_mod, newdata = test_df)
  rss_null  <- sum((test_df$phenotype - pred_null)^2)
  
  # 2. MIXPRS 
  mix_mod  <- lm(as.formula(paste("phenotype ~ MIXPRS +", paste(cov_choice, collapse="+"))), data = train_df) 
  pred_mix <- predict(mix_mod, newdata = test_df)
  rss_mix  <- sum((test_df$phenotype - pred_mix)^2)
  mix_r2_folds[i] <- 1 - (rss_mix / rss_null)
  
  # 3. OLS
  ols_formula <- as.formula(paste("phenotype ~", paste(c(cov_choice, indprs_col_names), collapse="+")))
  ind_mod  <- lm(ols_formula, data = train_df)
  pred_ind <- predict(ind_mod, newdata = test_df)
  rss_ind  <- sum((test_df$phenotype - pred_ind)^2)
  ind_r2_folds[i] <- 1 - (rss_ind / rss_null)
  
  # 4. Lasso
  cv_lasso <- cv.glmnet(x_train, y_train, alpha = 1, penalty.factor = p_factor, type.measure = "mse", nfolds = 4)
  pred_lasso <- predict(cv_lasso, s = "lambda.min", newx = x_test)
  rss_lasso  <- sum((y_test - as.vector(pred_lasso))^2)
  lasso_r2_folds[i] <- 1 - (rss_lasso / rss_null)
  
  # 5. Ridge
  cv_ridge <- cv.glmnet(x_train, y_train, alpha = 0, penalty.factor = p_factor, type.measure = "mse", nfolds = 4)
  pred_ridge <- predict(cv_ridge, s = "lambda.min", newx = x_test)
  rss_ridge  <- sum((y_test - as.vector(pred_ridge))^2)
  ridge_r2_folds[i] <- 1 - (rss_ridge / rss_null)
  
  # 6. Elastic Net
  alpha_seq <- seq(0.1, 0.9, 0.1)
  best_enet_cv_mse <- Inf
  best_enet_model <- NULL
  
  for (curr_alpha in alpha_seq) {
    cv_curr <- cv.glmnet(x_train, y_train, alpha = curr_alpha, penalty.factor = p_factor, type.measure = "mse", nfolds = 4)
    if (min(cv_curr$cvm) < best_enet_cv_mse) {
      best_enet_cv_mse <- min(cv_curr$cvm)
      best_enet_model <- cv_curr
    }
  }
  pred_enet <- predict(best_enet_model, s = "lambda.min", newx = x_test)
  rss_enet  <- sum((y_test - as.vector(pred_enet))^2)
  enet_r2_folds[i] <- 1 - (rss_enet / rss_null)
  
  # 7. NNLS
  y_tr_vec <- train_df$phenotype
  cov_tr_df <- train_df[, ..cov_choice]
  prs_tr_mat <- as.matrix(train_df[, ..indprs_col_names])
  prs_tr_mat <- apply(prs_tr_mat, 2, as.numeric)
  
  fit_y_cov <- lm(y_tr_vec ~ ., data = cov_tr_df)
  y_tr_res <- residuals(fit_y_cov)
  
  prs_tr_res <- matrix(NA, nrow = nrow(prs_tr_mat), ncol = ncol(prs_tr_mat))
  colnames(prs_tr_res) <- colnames(prs_tr_mat)
  for (j in 1:ncol(prs_tr_mat)) {
    temp_df <- cbind(cov_tr_df, prs = prs_tr_mat[, j])
    fit_prs_j <- lm(prs ~ ., data = temp_df)
    prs_tr_res[, j] <- residuals(fit_prs_j)
  }
  
  # NNLS on residuals (already standardized PRS, residuals should be comparable)
  nnls_fit <- nnls(as.matrix(prs_tr_res), y_tr_res)
  nnls_weights <- coef(nnls_fit)
  
  # Calculate scores using standardized PRS
  prs_tr_mat_original <- as.matrix(train_df[, ..indprs_col_names])
  prs_tr_mat_original <- apply(prs_tr_mat_original, 2, as.numeric)
  nnls_score_train <- prs_tr_mat_original %*% nnls_weights
  
  prs_te_mat_original <- as.matrix(test_df[, ..indprs_col_names])
  prs_te_mat_original <- apply(prs_te_mat_original, 2, as.numeric)
  nnls_score_test  <- prs_te_mat_original %*% nnls_weights
  
  train_df_nnls <- train_df
  train_df_nnls$NNLS_SCORE <- as.numeric(nnls_score_train)
  test_df_nnls <- test_df
  test_df_nnls$NNLS_SCORE <- as.numeric(nnls_score_test)
  
  final_nnls_mod <- lm(as.formula(paste("phenotype ~ NNLS_SCORE +", paste(cov_choice, collapse="+"))), 
                       data = train_df_nnls)
  pred_nnls <- predict(final_nnls_mod, newdata = test_df_nnls)
  
  rss_nnls <- sum((test_df$phenotype - pred_nnls)^2)
  nnls_r2_folds[i] <- 1 - (rss_nnls / rss_null)
}
      
      # All
      new_row <- data.frame(
        group = group,
        trait = trait,
        population = pop,
        n_samples = n_samples,
        n_prs_scores = length(indprs_col_names),
        
        # Means
        mixprs_r2_mean = mean(mix_r2_folds),
        indprs_r2_mean = mean(ind_r2_folds),
        lasso_r2_mean  = mean(lasso_r2_folds),
        ridge_r2_mean  = mean(ridge_r2_folds),
        enet_r2_mean   = mean(enet_r2_folds),
        nnls_r2_mean   = mean(nnls_r2_folds),
        
        # Folds - OLS
        indprs_r2_fold1 = ind_r2_folds[1],
        indprs_r2_fold2 = ind_r2_folds[2],
        indprs_r2_fold3 = ind_r2_folds[3],
        indprs_r2_fold4 = ind_r2_folds[4],
        
        # Folds - LassoPRS
        lasso_r2_fold1 = lasso_r2_folds[1],
        lasso_r2_fold2 = lasso_r2_folds[2],
        lasso_r2_fold3 = lasso_r2_folds[3],
        lasso_r2_fold4 = lasso_r2_folds[4],
        
        # Folds - RidgePRS
        ridge_r2_fold1 = ridge_r2_folds[1],
        ridge_r2_fold2 = ridge_r2_folds[2],
        ridge_r2_fold3 = ridge_r2_folds[3],
        ridge_r2_fold4 = ridge_r2_folds[4],
        
        # Folds - ElasticNet
        enet_r2_fold1 = enet_r2_folds[1],
        enet_r2_fold2 = enet_r2_folds[2],
        enet_r2_fold3 = enet_r2_folds[3],
        enet_r2_fold4 = enet_r2_folds[4],
        
        # Folds - NNLS
        nnls_r2_fold1 = nnls_r2_folds[1],
        nnls_r2_fold2 = nnls_r2_folds[2],
        nnls_r2_fold3 = nnls_r2_folds[3],
        nnls_r2_fold4 = nnls_r2_folds[4],
        
        # Folds - MIXPRS
        mixprs_r2_fold1 = mix_r2_folds[1],
        mixprs_r2_fold2 = mix_r2_folds[2],
        mixprs_r2_fold3 = mix_r2_folds[3],
        mixprs_r2_fold4 = mix_r2_folds[4],
        
        stringsAsFactors = FALSE
      )
      
      results_df <- rbind(results_df, new_row)
    }
  }
}

final_csv_path <- file.path("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/evaluation/IndPRS_Lasso_Ridge_Enet_NNLS_vs_MIXPRS_Evaluation.csv")
write.csv(results_df, final_csv_path, row.names = FALSE, quote = FALSE)