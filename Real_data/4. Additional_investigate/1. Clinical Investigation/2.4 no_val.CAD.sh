library(data.table)
library(stringr)
library(dplyr)
library(pROC)
library(ResourceSelection)

cov_choice = c("age_recruit","sex",paste0("PC",1:20))
total_covariates = fread("/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/cov_data/adjust.csv")
total_covariates = total_covariates %>% select(all_of(c("eid", cov_choice)))

pop = "EAS"
trait = "CAD"

## test phenotype
Trait_pheno <- fread(paste0("/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/pheno_data/",trait,"/",trait,"_",pop,".tsv"))
colnames(Trait_pheno) = c("eid","pheno")
Trait_pheno_id = Trait_pheno$eid

## PRS scores aligned with phenotype
Trait_ClinicalPRS = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/clinicalPRS/",trait,"/UKB_",trait,"_ClinicalPRS_prs_",pop,".sscore"))
Trait_ClinicalPRS = Trait_ClinicalPRS[match(Trait_pheno_id,Trait_ClinicalPRS$IID)]

Trait_SBayesRC_EUR = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/clinicalPRS/",trait,"/UKB_",trait,"_EUR_inter_SBayesRC_HM3_prs_",pop,".sscore"))
Trait_SBayesRC_EUR = Trait_SBayesRC_EUR[match(Trait_pheno_id,Trait_SBayesRC_EUR$IID)]

Trait_SBayesRC_EAS = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/clinicalPRS/",trait,"/UKB_",trait,"_EAS_inter_SBayesRC_HM3_prs_",pop,".sscore"))
Trait_SBayesRC_EAS = Trait_SBayesRC_EAS[match(Trait_pheno_id,Trait_SBayesRC_EAS$IID)]

Trait_MIXPRS = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/MIXPRS/UKB_",trait,"_MIXPRS_subsample_prune_prs_",pop,".sscore"))
Trait_MIXPRS = Trait_MIXPRS[match(Trait_pheno_id,Trait_MIXPRS$IID),]

Trait_MIXPRS_plus = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/clinicalPRS/",trait,"/UKB_",trait,"_MIXPRS_plus_prs_",pop,".sscore"))
Trait_MIXPRS_plus = Trait_MIXPRS_plus[match(Trait_pheno_id,Trait_MIXPRS_plus$IID)]

Trait_MIXPRS_plus_ClinicalPRS_SBayesRC = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/clinicalPRS/",trait,"/UKB_",trait,"_MIXPRS_plus_ClinicalPRS_SBayesRC_prs_",pop,".sscore"))
Trait_MIXPRS_plus_ClinicalPRS_SBayesRC = Trait_MIXPRS_plus_ClinicalPRS_SBayesRC[match(Trait_pheno_id,Trait_MIXPRS_plus_ClinicalPRS_SBayesRC$IID)]

## phenotype only
pheno = Trait_pheno[,-1]

method_names <- c("ClinicalPRS","SBayesRC_EUR","SBayesRC_EAS",
                  "MIXPRS","MIXPRS_plus","MIXPRS_plus_ClinicalPRS_SBayesRC")

test_data = data.table(
  ClinicalPRS = as.numeric(scale(Trait_ClinicalPRS$SCORE1_AVG)),
  SBayesRC_EUR  = as.numeric(scale(Trait_SBayesRC_EUR$SCORE1_AVG)),
  SBayesRC_EAS  = as.numeric(scale(Trait_SBayesRC_EAS$SCORE1_AVG)),
  MIXPRS = as.numeric(scale(Trait_MIXPRS$SCORE1_AVG)), 
  MIXPRS_plus = as.numeric(scale(Trait_MIXPRS_plus$SCORE1_AVG)), 
  MIXPRS_plus_ClinicalPRS_SBayesRC = as.numeric(scale(Trait_MIXPRS_plus_ClinicalPRS_SBayesRC$SCORE1_AVG))
)

cov_sub <- total_covariates[match(Trait_pheno_id, total_covariates$eid)]
PCnames <- paste0("PC", 1:20)

analysis_dt <- data.table(
  eid   = Trait_pheno_id,
  pheno = pheno$pheno
)
analysis_dt <- cbind(
  analysis_dt,
  cov_sub[, ..cov_choice],  # age_recruit, sex, PC1..PC20
  test_data                 # 6 PRS columns (scaled)
)

## 1) AUC (PRS only)
AUC_out <- data.table()

for (j in seq_along(method_names)) {
  Method <- method_names[j]
  data   <- pheno
  data$prs <- unlist(test_data[, ..j])

  glmfit      <- glm(pheno ~ prs, data = data, family = binomial(link = "logit"))
  glmfit_prob <- predict(glmfit, type = "response")
  glmfit_auc  <- roc(data$pheno, glmfit_prob, quiet = TRUE, plot = FALSE)$auc

  AUC_out <- rbind(
    AUC_out,
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
  AUC_out,
  paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/evaluation/clinical_",trait,"_",pop,"_AUC.csv"),
  quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE
)

## 2) OR per SD (adjusted)
OR_per_SD_out <- data.table()

for (Method in method_names) {
  dt <- copy(analysis_dt)
  dt <- dt[complete.cases(pheno, age_recruit, sex, dt[[Method]], dt[, ..PCnames])]
  
  dt[, pred := as.numeric(scale(get(Method)))]

  fml <- as.formula(
    paste("pheno ~ pred +", paste(cov_choice, collapse = " + "))
  )
  fit <- glm(fml, data = dt, family = binomial())

  co <- summary(fit)$coefficients["pred", ]

  OR_per_SD_out <- rbind(
    OR_per_SD_out,
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
  OR_per_SD_out,
  paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/evaluation/clinical_",trait,"_",pop,"_OR_perSD.csv"),
  quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE
)

## 3) OR by 5 quantile groups (each vs bottom 20%)
make_quantile5 <- function(x) {
  d <- dplyr::ntile(x, 5)  # 1..5
  factor(
    d, levels = 1:5,
    labels = c("<20%","20-40%","40-60%","60-80%",">=80%")
  )
}

quantile5_or_out <- data.table()

for (Method in method_names) {
  dt <- copy(analysis_dt)
  dt <- dt[complete.cases(pheno, age_recruit, sex, dt[[Method]], dt[, ..PCnames])]

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

  quantile5_or_out <- rbind(quantile5_or_out, out_rows, fill = TRUE)
}

grp_levels <- c("<20%","20-40%","40-60%","60-80%",">=80%")
quantile5_or_out[, Group := factor(Group, levels = grp_levels)]

write.table(
  quantile5_or_out,
  paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/evaluation/clinical_",trait,"_",pop,"_OR_by_quantile5_vs_bottom20.csv"),
  quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE
)

## 4) Calibration (HL test + quintile summaries)
calibration_pval <- data.table()
calibration_out  <- data.table()

for (Method in method_names) {
  tmp <- copy(analysis_dt)
  tmp <- tmp[complete.cases(pheno, age_recruit, sex, tmp[[Method]], tmp[, ..PCnames])]

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
    calibration_pval <- rbind(
      calibration_pval,
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
    calibration_pval <- rbind(
      calibration_pval,
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

  calibration_out <- rbind(calibration_out, dec_tbl, fill = TRUE)
}

write.table(
  calibration_pval,
  paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/evaluation/clinical_",trait,"_",pop,"_calibration_pval_pred.csv"),
  quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE
)
write.table(
  calibration_out,
  paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/evaluation/clinical_",trait,"_",pop,"_calibration_forPlot_pred.csv"),
  quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE
)

## 5) Prevalence per PRS quintile
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
    ci$prevalence <- ci$prevalence[order(ci$quantile_group)]
  } else {
    ci$quantile_group <- ci$Group.1
  }

  data.table(
    quantile_group = ci$quantile_group,
    prevalence     = ci$prevalence
  )
}

prev_out <- data.table()

for (Method in method_names) {
  target <- analysis_dt$pheno
  prs    <- analysis_dt[[Method]]
  idx    <- complete.cases(target, prs)

  add <- getPrevByQuantile(target[idx], prs[idx])
  add <- cbind(
    pop    = pop,
    trait  = trait,
    Method = Method,
    add
  )

  prev_out <- rbind(prev_out, add, fill = TRUE)
}

write.table(
  prev_out,
  paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/evaluation/clinical_",trait,"_",pop,"_prevalence_per_quantile.csv"),
  quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE
)