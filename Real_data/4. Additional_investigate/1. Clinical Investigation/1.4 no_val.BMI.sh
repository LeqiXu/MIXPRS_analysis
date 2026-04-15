library(data.table)
library(stringr)
library(dplyr)
library(pROC)
library(ResourceSelection)

cov_choice <- c("age_recruit","sex",paste0("PC",1:20))
total_covariates <- fread("/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/cov_data/adjust.csv")
total_covariates <- total_covariates %>% dplyr::select(all_of(c("eid", cov_choice)))

trait <- "BMI"
pops  <- c("EAS","AFR")

# Methods for BMI
method_names <- c("ClinicalPRS_META","SBayesRC_EUR","SBayesRC_EAS","SBayesRC_AFR",
                  "MIXPRS","MIXPRS_plus","MIXPRS_plus_ClinicalPRS_SBayesRC")
PCnames <- paste0("PC", 1:20)

make_quantile5 <- function(x) {
  d <- dplyr::ntile(x, 5)  # 1..5
  factor(
    d, levels = 1:5,
    labels = c("<20%","20-40%","40-60%","60-80%",">=80%")
  )
}

getPrevByQuantile <- function(y, x) {
  qbreaks <- quantile(x, probs = c(seq(0, 0.8, 0.2), 1),
                      na.rm = TRUE, names = FALSE)
  qbreaks <- unique(qbreaks)
  quantile_group <- cut(x, breaks = qbreaks, include.lowest = TRUE)

  ci <- aggregate(y, list(quantile_group), mean, na.rm = TRUE)
  colnames(ci)[2] <- "prevalence"

  intended_labels <- c("<20%","20-40%","40-60%","60-80%",">=80%")
  if (nrow(ci) == length(intended_labels)) {
    ci$quantile_group <- factor(intended_labels, levels = intended_labels)
    ci$prevalence     <- ci$prevalence[order(ci$quantile_group)]
  } else {
    ci$quantile_group <- ci$Group.1
  }

  data.table(
    quantile_group = ci$quantile_group,
    prevalence     = ci$prevalence
  )
}

out_dir <- "/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/evaluation/"

for (pop in pops) {
  cat("Processing BMI for", pop, "...\n")

  ## 1) Phenotype: continuous BMI + obesity
  Trait_pheno <- fread(paste0("/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/pheno_data/",
                              trait,"/",trait,"_",pop,".tsv"))
  colnames(Trait_pheno) <- c("eid","BMI")

  # Define obesity by ancestry
  if (pop %in% c("EAS","SAS")) {
    Trait_pheno[, obesity := ifelse(BMI >= 27.5, 1, 0)]
  } else {
    Trait_pheno[, obesity := ifelse(BMI >= 30, 1, 0)]
  }

  Trait_pheno_id <- Trait_pheno$eid

  # Continuous BMI (scaled) for R²
  pheno_cont <- as.numeric(scale(Trait_pheno$BMI))
  # Binary obesity for logistic metrics
  pheno_bin  <- Trait_pheno$obesity

  ## 2) PRS scores aligned with phenotype
  Trait_ClinicalPRS_META <- fread(
    paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/clinicalPRS/",
           trait,"/UKB_",trait,"_ClinicalPRS_META_prs_",pop,".sscore")
  )
  Trait_ClinicalPRS_META <-
    Trait_ClinicalPRS_META[match(Trait_pheno_id, Trait_ClinicalPRS_META$IID)]

  Trait_SBayesRC_EUR <- fread(
    paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/clinicalPRS/",
           trait,"/UKB_",trait,"_EUR_inter_SBayesRC_HM3_prs_",pop,".sscore")
  )
  Trait_SBayesRC_EUR <-
    Trait_SBayesRC_EUR[match(Trait_pheno_id, Trait_SBayesRC_EUR$IID)]

  Trait_SBayesRC_EAS <- fread(
    paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/clinicalPRS/",
           trait,"/UKB_",trait,"_EAS_inter_SBayesRC_HM3_prs_",pop,".sscore")
  )
  Trait_SBayesRC_EAS <-
    Trait_SBayesRC_EAS[match(Trait_pheno_id, Trait_SBayesRC_EAS$IID)]

  Trait_SBayesRC_AFR <- fread(
    paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/clinicalPRS/",
           trait,"/UKB_",trait,"_AFR_inter_SBayesRC_HM3_prs_",pop,".sscore")
  )
  Trait_SBayesRC_AFR <-
    Trait_SBayesRC_AFR[match(Trait_pheno_id, Trait_SBayesRC_AFR$IID)]

  Trait_MIXPRS <- fread(
    paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/MIXPRS/UKB_",
           trait,"_MIXPRS_subsample_prune_prs_",pop,".sscore")
  )
  Trait_MIXPRS <- Trait_MIXPRS[match(Trait_pheno_id, Trait_MIXPRS$IID)]

  Trait_MIXPRS_plus <- fread(
    paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/clinicalPRS/",
           trait,"/UKB_",trait,"_MIXPRS_plus_prs_",pop,".sscore")
  )
  Trait_MIXPRS_plus <-
    Trait_MIXPRS_plus[match(Trait_pheno_id, Trait_MIXPRS_plus$IID)]

  Trait_MIXPRS_plus_ClinicalPRS_SBayesRC <- fread(
    paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/clinicalPRS/",
           trait,"/UKB_",trait,"_MIXPRS_plus_ClinicalPRS_SBayesRC_prs_",pop,".sscore")
  )
  Trait_MIXPRS_plus_ClinicalPRS_SBayesRC <-
    Trait_MIXPRS_plus_ClinicalPRS_SBayesRC[
      match(Trait_pheno_id, Trait_MIXPRS_plus_ClinicalPRS_SBayesRC$IID)
    ]

  test_data <- data.table(
    ClinicalPRS_META = as.numeric(scale(Trait_ClinicalPRS_META$SCORE1_AVG)),
    SBayesRC_EUR = as.numeric(scale(Trait_SBayesRC_EUR$SCORE1_AVG)),
    SBayesRC_EAS = as.numeric(scale(Trait_SBayesRC_EAS$SCORE1_AVG)),
    SBayesRC_AFR = as.numeric(scale(Trait_SBayesRC_AFR$SCORE1_AVG)),
    MIXPRS = as.numeric(scale(Trait_MIXPRS$SCORE1_AVG)),
    MIXPRS_plus = as.numeric(scale(Trait_MIXPRS_plus$SCORE1_AVG)),
    MIXPRS_plus_ClinicalPRS_SBayesRC = as.numeric(scale(Trait_MIXPRS_plus_ClinicalPRS_SBayesRC$SCORE1_AVG))
  )

  ## 3) Covariates aligned
  cov_sub <- total_covariates[match(Trait_pheno_id, total_covariates$eid)]

  ## 4) Analysis data.tables
  # continuous BMI
  analysis_dt_cont <- data.table(
    eid   = Trait_pheno_id,
    pheno = pheno_cont
  )
  analysis_dt_cont <- cbind(
    analysis_dt_cont,
    cov_sub[, ..cov_choice],
    test_data
  )

  # obesity (binary)
  analysis_dt_bin <- data.table(
    eid   = Trait_pheno_id,
    pheno = pheno_bin
  )
  analysis_dt_bin <- cbind(
    analysis_dt_bin,
    cov_sub[, ..cov_choice],
    test_data
  )

  ## (1) R² for continuous BMI (incremental over covariates)
  pheno_covariates <- data.frame(
    pheno = analysis_dt_cont$pheno,
    analysis_dt_cont[, ..cov_choice]
  )

  linear_null <- lm(pheno ~ ., data = pheno_covariates)
  linear_null_res2 <- sum(residuals(linear_null)^2)

  R2_out_pop <- data.table()

  for (j in seq_along(method_names)) {
    Method <- method_names[j]
    data_lin <- pheno_covariates
    data_lin$prs <- test_data[[j]]

    fit_lin <- lm(pheno ~ . + prs, data = data_lin)
    res2    <- sum(residuals(fit_lin)^2)
    R2_inc  <- 1 - res2 / linear_null_res2

    R2_out_pop <- rbind(
      R2_out_pop,
      data.table(
        pop    = pop,
        trait  = trait,
        Method = Method,
        R2     = R2_inc
      ),
      fill = TRUE
    )
  }

  write.table(
    R2_out_pop,
    paste0(out_dir,"clinical_",trait,"_",pop,"_r2.csv"),
    quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE
  )

  ## (2) Binary metric for obesity
  AUC_out_pop          <- data.table()
  OR_per_SD_out_pop    <- data.table()
  quantile5_or_out_pop <- data.table()
  calibration_pval_pop <- data.table()
  calibration_out_pop  <- data.table()
  prev_out_pop         <- data.table()

  grp_levels <- c("<20%","20-40%","40-60%","60-80%",">=80%")

  ## (2) AUC (PRS only) for obesity
  pheno_bin_df <- data.frame(pheno = pheno_bin)

  for (j in seq_along(method_names)) {
    Method <- method_names[j]
    data_auc <- pheno_bin_df
    data_auc$prs <- unlist(test_data[, ..j])

    glmfit      <- glm(pheno ~ prs, data = data_auc,
                       family = binomial(link = "logit"))
    glmfit_prob <- predict(glmfit, type = "response")
    glmfit_auc  <- roc(data_auc$pheno, glmfit_prob,
                       quiet = TRUE, plot = FALSE)$auc

    AUC_out_pop <- rbind(
      AUC_out_pop,
      data.table(
        pop    = pop,
        trait  = trait,
        Method = Method,
        AUC    = as.numeric(glmfit_auc)
      ),
      fill = TRUE
    )
  }

  write.table(
    AUC_out_pop,
    paste0(out_dir,"clinical_",trait,"_",pop,"_AUC.csv"),
    quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE
  )

  ## (3) OR per SD (PRS + covariates) for obesity
  for (Method in method_names) {
    dt <- copy(analysis_dt_bin)
    dt <- dt[complete.cases(
      pheno,
      age_recruit,
      sex,
      dt[, ..PCnames],
      dt[[Method]]
    )]

    dt[, pred := as.numeric(scale(get(Method)))]

    fml <- as.formula(
      paste("pheno ~ pred +", paste(cov_choice, collapse = " + "))
    )
    fit <- glm(fml, data = dt, family = binomial())

    co <- summary(fit)$coefficients["pred", ]

    OR_per_SD_out_pop <- rbind(
      OR_per_SD_out_pop,
      data.table(
        pop    = pop,
        trait  = trait,
        Method = Method,
        OR     = exp(co["Estimate"]),
        CI_low = exp(co["Estimate"] - 1.96 * co["Std. Error"]),
        CI_high= exp(co["Estimate"] + 1.96 * co["Std. Error"]),
        pval   = co["Pr(>|z|)"]
      ),
      fill = TRUE
    )
  }

  write.table(
    OR_per_SD_out_pop,
    paste0(out_dir,"clinical_",trait,"_",pop,"_OR_perSD.csv"),
    quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE
  )

  ## (4) OR by 5 PRS quantile groups (vs bottom 20%) for obesity
  for (Method in method_names) {
    dt <- copy(analysis_dt_bin)
    dt <- dt[complete.cases(
      pheno,
      age_recruit,
      sex,
      dt[, ..PCnames],
      dt[[Method]]
    )]

    dt[, group := make_quantile5(get(Method))]
    dt[, group := stats::relevel(group, ref = "<20%")]

    fml <- as.formula(
      paste("pheno ~ group +", paste(cov_choice, collapse = " + "))
    )
    fit <- glm(fml, data = dt, family = binomial())
    sm  <- summary(fit)$coefficients

    by_grp <- dt[, .N, by = group][order(group)]

    out_rows <- data.table(
      pop    = pop,
      trait  = trait,
      Method = Method,
      Group  = as.character(by_grp$group),
      n      = by_grp$N,
      OR     = NA_real_,
      CI_low = NA_real_,
      CI_high= NA_real_,
      pval   = NA_real_
    )

    coef_names <- rownames(sm)
    grp_rows   <- coef_names[grepl("^group", coef_names)]
    if (length(grp_rows) > 0) {
      tmp <- data.table(
        coef_name = grp_rows,
        Estimate  = sm[grp_rows, "Estimate"],
        SE        = sm[grp_rows, "Std. Error"],
        pval      = sm[grp_rows, "Pr(>|z|)"]
      )
      tmp[, OR      := exp(Estimate)]
      tmp[, CI_low  := exp(Estimate - 1.96 * SE)]
      tmp[, CI_high := exp(Estimate + 1.96 * SE)]
      tmp[, Group   := sub("^group", "", coef_name)]

      out_rows[Group != "<20%", `:=`(
        OR     = tmp$OR[match(Group, tmp$Group)],
        CI_low = tmp$CI_low[match(Group, tmp$Group)],
        CI_high= tmp$CI_high[match(Group, tmp$Group)],
        pval   = tmp$pval[match(Group, tmp$Group)]
      )]
    }

    out_rows[Group == "<20%", `:=`(OR = 1, CI_low = 1, CI_high = 1, pval = NA_real_)]
    out_rows[, Group := factor(Group, levels = grp_levels)]

    quantile5_or_out_pop <- rbind(quantile5_or_out_pop, out_rows, fill = TRUE)
  }

  write.table(
    quantile5_or_out_pop,
    paste0(out_dir,"clinical_",trait,"_",pop,"_OR_by_quantile5_vs_bottom20.csv"),
    quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE
  )

  ## (5) Calibration: HL test + quintile summaries for obesity
  for (Method in method_names) {
    tmp <- copy(analysis_dt_bin)
    tmp <- tmp[complete.cases(
      pheno,
      age_recruit,
      sex,
      tmp[, ..PCnames],
      tmp[[Method]]
    )]

    fml <- as.formula(
      paste("pheno ~ scale(", Method, ") + ",
            paste(cov_choice, collapse = " + "), sep = "")
    )
    fit <- glm(fml, data = tmp, family = binomial())
    tmp$predicted <- predict(fit, type = "response")

    hl_res <- tryCatch(
      hoslem.test(tmp$pheno, tmp$predicted, g = 5),
      error = function(e) NULL
    )

    if (!is.null(hl_res)) {
      calibration_pval_pop <- rbind(
        calibration_pval_pop,
        data.table(
          pop    = pop,
          trait  = trait,
          Method = Method,
          HL_stat= unname(hl_res$statistic["X-squared"]),
          HL_df  = unname(hl_res$parameter),
          HL_pval= unname(hl_res$p.value),
          g      = 5L
        ),
        fill = TRUE
      )
    } else {
      calibration_pval_pop <- rbind(
        calibration_pval_pop,
        data.table(
          pop    = pop,
          trait  = trait,
          Method = Method,
          HL_stat= NA_real_,
          HL_df  = NA_real_,
          HL_pval= NA_real_,
          g      = 5L
        ),
        fill = TRUE
      )
    }

    dec_tbl <- data.table(
      pop    = pop,
      trait  = trait,
      Method = Method,
      group  = dplyr::ntile(tmp$predicted, 5),
      pred   = tmp$predicted,
      obs    = tmp$pheno
    )[
      , .(mean_pred = mean(pred, na.rm = TRUE),
          mean_obs  = mean(obs,  na.rm = TRUE),
          n         = .N),
      by = .(pop, trait, Method, group)
    ][order(pop, trait, Method, group)]

    calibration_out_pop <- rbind(calibration_out_pop, dec_tbl, fill = TRUE)
  }

  write.table(
    calibration_pval_pop,
    paste0(out_dir,"clinical_",trait,"_",pop,"_calibration_pval_pred.csv"),
    quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE
  )

  write.table(
    calibration_out_pop,
    paste0(out_dir,"clinical_",trait,"_",pop,"_calibration_forPlot_pred.csv"),
    quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE
  )

  ## (6) Prevalence per PRS quintile for obesity
  for (Method in method_names) {
    target <- analysis_dt_bin$pheno
    prs    <- analysis_dt_bin[[Method]]
    idx    <- complete.cases(target, prs)

    add <- getPrevByQuantile(target[idx], prs[idx])
    add <- cbind(
      pop    = pop,
      trait  = trait,
      Method = Method,
      add
    )

    prev_out_pop <- rbind(prev_out_pop, add, fill = TRUE)
  }

  write.table(
    prev_out_pop,
    paste0(out_dir,"clinical_",trait,"_",pop,"_prevalence_per_quantile.csv"),
    quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE
  )
}

