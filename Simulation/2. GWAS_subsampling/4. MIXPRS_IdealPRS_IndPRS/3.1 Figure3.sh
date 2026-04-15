# ---- Setup ----
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(ggpubr)

setwd("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/plot")

my_theme <- theme(
  # shrink axis titles
  axis.title.x = element_text(size = 8, face = "bold"),
  axis.title.y = element_text(size = 8, face = "bold"),
  # shrink tick labels
  axis.text.x  = element_text(size = 5, face = "bold"),
  axis.text.y  = element_text(size = 6, face = "bold"),
  # leave all other text at your preferred size
  plot.title   = element_text(size = 10, face = "bold"),
  strip.text   = element_text(size = 10, face = "bold"),
  panel.grid.major.x = element_blank(),
  panel.grid.minor.x = element_blank(),
  panel.grid.major.y = element_blank(),
  panel.grid.minor.y = element_blank(),
  legend.position  = "top",
  legend.title     = element_text(size = 9, face = "bold"),     # title size
  legend.text      = element_text(size = 7),
  plot.margin = unit(c(0.1, 0.1, 0.5, 0.1), "cm")
)

compare_method_color <- c("IndPRS"  = "#0072B2", "IdealPRS" = "#E69F00", "MIXPRS"  = "#6A0DAD")
ind_method_color <- c("Ridge" = "#E69F00", "Lasso" = "#D55E00", "Elastic Net"= "#0072B2","NNLS" = "#6A0DAD")

## load data
prs_table <- fread("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/evaluation/sim_PRS_real_MIXPRS_IdealPRS_IndPRS_r2.csv")
prs_long_table <- melt(prs_table, id.vars = c("n", "pop", "p", "rhog", "sample1", "sample2"),
                   variable.name = "method", value.name = "r2")
prs_long_table$pop <- factor(prs_long_table$pop, levels = c("EUR","EAS","AFR","SAS","AMR"))
prs_long_table$p <- factor(prs_long_table$p,
                       levels = c(0.1, 0.01, 0.001, 5e-04),
                       labels = c("p = 0.1", "p = 0.01", "p = 0.001", "p = 5 × 10⁻⁴"))
prs_summary <- prs_long_table %>%
  group_by(pop, p, method) %>%
  summarise(mean_r2 = mean(r2), se_r2 = sd(r2)/sqrt(n()), .groups = "drop")

## p1: MIXPRS, IdealPRS, IndPRS
compare_summary = prs_summary[which(prs_summary$method %in% c("IndPRS_nnls","IdealPRS_nnls","MIXPRS_nnls")),]
compare_summary$method <- factor(compare_summary$method, levels = c("IndPRS_nnls", "IdealPRS_nnls", "MIXPRS_nnls"),
                                                         labels = c("IndPRS", "IdealPRS", "MIXPRS"))
compare_summary$method <- factor(compare_summary$method, levels = c("IndPRS", "IdealPRS", "MIXPRS"))
setorder(compare_summary, method)

p1 <- ggplot(compare_summary, aes(x = pop, y = mean_r2, color = method, group = method)) +
  geom_point(size = 1) +
  facet_wrap(~ p, ncol = 5) +
  labs(title = "",
       x = "",
       y = "R²",
       color = "Method  ") +
  theme_classic() +
  my_theme +
  scale_color_manual(values = compare_method_color) +
  scale_y_continuous(limits = c(NA, 0.41))

## p2: IndPRS_nnls, IndPRS_ridge, IndPRS_lasso, IndPRS_elasticnet
ind_summary = prs_summary[which(prs_summary$method %in% c("IndPRS_ridge","IndPRS_lasso","IndPRS_elasticnet","IndPRS_nnls")),]
ind_summary$method = factor(ind_summary$method, levels = c("IndPRS_ridge","IndPRS_lasso","IndPRS_elasticnet","IndPRS_nnls"),
                                                labels = c("Ridge","Lasso","Elastic Net","NNLS"))
setorder(ind_summary, method)

p2 <- ggplot(ind_summary, aes(x = pop, y = mean_r2, color = method, group = method)) +
  geom_point(size = 1) +
  facet_wrap(~ p, ncol = 5) +
  labs(title = "",
       x = "",
       y = "R²",
       color = "Strategy  ") +
  theme_classic() +
  my_theme +
  scale_color_manual(values = ind_method_color) +
  scale_y_continuous(limits = c(NA, 0.41))

## combine plot
combined_plot <- ggarrange(p1, p2, ncol = 1, nrow = 2, labels = c("a","b"), heights = c(0.5,0.5))
print(combined_plot)

# Save figure (change path/size as needed)
ggsave(filename = "Figure3.pdf",
       plot = combined_plot,
       width = 6.5, 
       height = 8,
       units = "in",
       dpi = 300,
       device   = cairo_pdf) 