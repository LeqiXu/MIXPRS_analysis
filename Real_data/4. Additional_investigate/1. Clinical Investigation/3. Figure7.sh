library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)

my_theme <- theme(
  # shrink axis titles
  axis.title.x = element_text(size = 8, face = "bold"),
  axis.title.y = element_text(size = 8, face = "bold"),
  # shrink tick labels
  axis.text.x  = element_text(size = 5, face = "bold"),
  axis.text.y  = element_text(size = 6, face = "bold"),
  # leave all other text at your preferred size
  plot.title   = element_text(size = 8, face = "bold", hjust = 0.5),
  strip.text   = element_text(size = 8, face = "bold"),
  panel.grid.major.x = element_blank(),
  panel.grid.minor.x = element_blank(),
  panel.grid.major.y = element_blank(),
  panel.grid.minor.y = element_blank(),
  legend.position  = "bottom",
  legend.justification = "center",
  legend.box.just = "center",
  legend.title     = element_text(size = 6, face = "bold"),     # title size
  legend.text      = element_text(size = 6),
  plot.margin = unit(c(0.1, 0.1, 0.5, 0.1), "cm")
)

setwd("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/plot")

eval_dir <- "/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/evaluation/"

make_trait_plot <- function(trait, pop, include_R2 = FALSE, eval_dir = NULL) {

  message("Processing ", trait, " - ", pop, " (include_R2 = ", include_R2, ") ...")

  ## ---------- Read core evaluation files ----------
  AUC_dt   <- fread(paste0(eval_dir, "clinical_", trait, "_", pop, "_AUC.csv"))
  ORsd_dt  <- fread(paste0(eval_dir, "clinical_", trait, "_", pop, "_OR_perSD.csv"))
  ORq_dt   <- fread(paste0(eval_dir, "clinical_", trait, "_", pop, "_OR_by_quantile5_vs_bottom20.csv"))
  CAL_dt   <- fread(paste0(eval_dir, "clinical_", trait, "_", pop, "_calibration_pval_pred.csv"))

  ## R2 (only if include_R2 = TRUE)
  R2_dt <- NULL
  if (include_R2) {
    r2_file <- paste0(eval_dir, "clinical_", trait, "_", pop, "_r2.csv")
    if (!file.exists(r2_file)) {
      stop("include_R2 = TRUE but R2 file not found: ", r2_file)
    }
    R2_dt <- fread(r2_file)
  }

  ## ---------- Normalize Method names ----------
  norm_method <- function(dt) {
    if (!is.null(dt) && "Method" %in% names(dt)) {
      dt[Method == "ClinicalPRS_META", Method := "ClinicalPRS"]
    }
    dt
  }

  AUC_dt  <- norm_method(AUC_dt)
  ORsd_dt <- norm_method(ORsd_dt)
  ORq_dt  <- norm_method(ORq_dt)
  CAL_dt  <- norm_method(CAL_dt)
  R2_dt   <- norm_method(R2_dt)

  ## ---------- Summaries ----------
  AUC_sum   <- AUC_dt[,  .(Method, AUC)]
  ORsd_sum  <- ORsd_dt[, .(Method, OR_perSD = OR)]
  ORtop_sum <- ORq_dt[Group == ">=80%", .(Method, OR_top20 = OR)]
  CAL_sum   <- CAL_dt[, .(Method, HL_pval)]

  if (!is.null(R2_dt)) {
    R2_sum <- R2_dt[, .(Method, R2)]
  }

  ## ---------- Merge ----------
  summary_dt <- AUC_sum %>%
    inner_join(ORsd_sum,  by = "Method") %>%
    inner_join(ORtop_sum, by = "Method") %>%
    inner_join(CAL_sum,   by = "Method")

  if (!is.null(R2_dt)) {
    summary_dt <- summary_dt %>% inner_join(R2_sum, by = "Method")
  }

  summary_dt <- as.data.table(summary_dt)

  ## ---------- Collapse SBayesRC_* → SBayesRC (max across metrics) ----------
  sbrc_list <- grep("^SBayesRC", summary_dt$Method, value = TRUE)

  if (length(sbrc_list) > 1) {
    metric_cols <- setdiff(names(summary_dt), "Method")
    sbrc_dt <- summary_dt[Method %in% sbrc_list]
    max_vals <- sbrc_dt[, lapply(.SD, max, na.rm = TRUE), .SDcols = metric_cols]
    SBayesRC_combined <- data.table(Method = "SBayesRC", max_vals)

    summary_dt <- summary_dt[!Method %in% sbrc_list]
    summary_dt <- rbind(summary_dt, SBayesRC_combined, fill = TRUE)
  }

  ## ---------- Keep only main methods ----------
  method_list     <- c("MIXPRS","SBayesRC","ClinicalPRS","MIXPRS_plus_ClinicalPRS_SBayesRC")
  method_namelist <- c("MIXPRS","SBayesRC","ClinicalPRS","MIXPRS+")

  summary_dt <- summary_dt[Method %in% method_list]
  summary_dt$Method <- factor(summary_dt$Method, levels = method_list, labels = method_namelist)

  ## ---------- Ranking ----------
  if (include_R2) {
    summary_dt[, R2_rank := rank(-R2)]
  }

  summary_dt[, AUC_rank    := rank(-AUC)]
  summary_dt[, OR_SD_rank  := rank(-OR_perSD)]
  summary_dt[, OR_top_rank := rank(-OR_top20)]
  summary_dt[, HL_rank     := rank(-HL_pval)]

  if (include_R2) {
    summary_dt[, overall_rank := (R2_rank + AUC_rank + OR_SD_rank + OR_top_rank + HL_rank) / 5]
    metric_cols_for_plot <- c("overall_rank","R2","AUC","OR_perSD","OR_top20","HL_pval")
    metric_labels <- c("Average rank","R²","AUC",
                       "OR per SD","OR (>=80% vs <20%)","HL p-value")
  } else {
    summary_dt[, overall_rank := (AUC_rank + OR_SD_rank + OR_top_rank + HL_rank) / 4]
    metric_cols_for_plot <- c("overall_rank","AUC","OR_perSD","OR_top20","HL_pval")
    metric_labels <- c("Average rank","AUC",
                       "OR per SD","OR (>=80% vs <20%)","HL p-value")
  }

  ## ---------- Long format ----------
  plot_dt <- melt(
    summary_dt,
    id.vars = c("Method","overall_rank"),
    measure.vars = metric_cols_for_plot,
    variable.name = "Metric",
    value.name    = "Value"
  )

  plot_dt[, Metric := factor(Metric, levels = metric_cols_for_plot, labels = metric_labels)]

  plot_dt[, is_best_metric := {
    if (Metric[1] == "Average rank") Value == min(Value, na.rm = TRUE)
    else Value == max(Value, na.rm = TRUE)
  }, by = Metric]

  ## ---------- Overall plot ----------
  overall_dt <- plot_dt[Metric == "Average rank"]
  overall_dt[, Method := factor(Method, levels = overall_dt[order(Value)]$Method)]

  gg_overall <- ggplot(overall_dt, aes(Method, Value)) +
    geom_segment(aes(xend = Method, y = 0, yend = Value, color = is_best_metric)) +
    geom_point(aes(color = is_best_metric), size = 1.5) +
    scale_color_manual(
      values = c("FALSE" = "grey50", "TRUE" = "#08519c"),
      labels = c("Not best", "Best"),
      name   = "Best average rank (lower = better)"
    ) +
    coord_flip() +
    guides(color = guide_legend(
      direction      = "horizontal",
      title.position = "top",
      title.hjust    = 0.5,
      nrow           = 1
    )) +
    labs(x = "",
         y = "Average rank",
         title = paste0("UKBB ", trait," ",pop)) +
    theme_classic() + my_theme

  ## ---------- Heatmap ----------
  metric_dt <- plot_dt[Metric != "Average rank"]
  overall_order <- overall_dt[order(Value), Method]
  metric_dt[, Method := factor(Method, levels = overall_order)]

  metric_dt[, scaled_value := {
    r <- range(Value, na.rm = TRUE)
    if (diff(r) == 0) 0.5 else (Value - r[1]) / diff(r)
  }, by = Metric]

  gg_metric <- ggplot(metric_dt, aes(Metric, Method)) +
    geom_tile(aes(fill = scaled_value), color = "white") +
    geom_text(aes(label = round(Value, 3),
                  fontface = ifelse(is_best_metric, "bold","plain")),
                  size = 2) +
    scale_fill_gradientn(
      colours = c("#f7fbff", "#d7e6f6", "#aecfe6", "#7fb6d8", "#4f97c8"),
      name    = "Relative metric performance (higher = better)"
    ) +
    guides(fill = guide_colorbar(
      title.position = "top",
      title.hjust    = 0.5, 
      label.position = "bottom",
      direction      = "horizontal",
      barheight      = unit(0.2, "cm")
    )) +
    labs(x = "",
         y = "",
         title = paste0("UKBB ", trait," ",pop)) +
    theme_classic() + my_theme

  list(overall = gg_overall,
       metric  = gg_metric)
} 

## ===== Generate your three plots =====

# CAD: no R2
CAD_EAS   <- make_trait_plot(trait = "CAD", pop = "EAS", include_R2 = FALSE, eval_dir = eval_dir)

CAD_EAS_overall_plot  <- CAD_EAS$overall
CAD_EAS_metric_plot   <- CAD_EAS$metric

# BMI: with R2
BMI_EAS   <- make_trait_plot(trait = "BMI", pop = "EAS", include_R2 = TRUE,  eval_dir = eval_dir)
BMI_AFR   <- make_trait_plot(trait = "BMI", pop = "AFR", include_R2 = TRUE,  eval_dir = eval_dir)

BMI_EAS_overall_plot  <- BMI_EAS$overall
BMI_EAS_metric_plot   <- BMI_EAS$metric

BMI_AFR_overall_plot  <- BMI_AFR$overall
BMI_AFR_metric_plot   <- BMI_AFR$metric

# overall plots
overall_plot <- ggarrange(
  BMI_EAS_overall_plot,
  BMI_AFR_overall_plot,
  CAD_EAS_overall_plot,
  nrow = 3,
  common.legend = TRUE,
  legend = "bottom"
)

# metric heatmaps
metric_plot <- ggarrange(
  BMI_EAS_metric_plot,
  BMI_AFR_metric_plot,
  CAD_EAS_metric_plot,
  nrow = 3,
  common.legend = TRUE,
  legend = "bottom"
)

# Combine rows
combined_plot <- ggarrange(
  overall_plot,
  metric_plot,
  ncol = 2,
  widths = c(0.35, 0.65),
  labels = c("a", "b")
)

print(combined_plot)

# Save figure (change path/size as needed)
ggsave(filename = "Figure7.pdf",
       plot = combined_plot,
       width = 6.5, 
       height = 8,
       units = "in",
       dpi = 300,
       device   = cairo_pdf) 