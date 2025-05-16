## MIXPRS comparison
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(gridExtra)
library(cowplot)
library(grid)
library(scales)
library(readr)
library(ggsci)

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
  panel.grid.major.y = element_line(color = "grey", size = 0.5),
  panel.grid.minor.y = element_line(color = "lightgrey", size = 0.25),
  legend.position  = "top",
  legend.title     = element_text(size = 9, face = "bold"),     # title size
  legend.text      = element_text(size = 7),
  plot.margin = unit(c(0.1, 0.1, 0.5, 0.1), "cm")
)

### Step1: Plot the correlation
library(data.table)
library(ggplot2)
library(dplyr)

corr_table = fread("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/residual/simulation_residual_corr.txt")
corr_table$pop = factor(corr_table$pop, levels = c("EUR","EAS","AFR","SAS","AMR"))
corr_table$type <- dplyr::recode(
  corr_table$type,
  "full_snplist"  = "RefLD_Full",
  "prune_snplist" = "RefLD_Prune",
  "prune_snplist_ind_approx"  = "Identity_Prune"
)
corr_table$type = factor(corr_table$type, levels = c("RefLD_Full","RefLD_Prune","Identity_Prune"))
corr_table$p <- factor(corr_table$p,
                       levels = c(0.1, 0.01, 0.001, 5e-04),
                       labels = c("p = 0.1", "p = 0.01", "p = 0.001", "p = 5 × 10⁻⁴"))

summary_corr_table<- corr_table %>%
  group_by(pop, p, type) %>%
  summarise(mean_corr = mean(corr), se_corr = sd(corr)/sqrt(n()), .groups = "drop")

p1 <- ggplot(summary_corr_table, aes(x = pop, y =mean_corr, color = type)) +
  geom_point(size = 1) +
  geom_line(aes(group = type)) + 
  geom_hline(yintercept = 0, linetype = "solid", color = alpha("darkred", alpha=0.5), size = 0.6) +
  facet_wrap(~ p, ncol = 4) +
  labs(x = NULL, y = "Residual Correlation", color = "Method") +
  theme_classic() +
  my_theme +
  scale_y_continuous(limits = c(-1, 1))

### Step2: Plot the performance
prs_table = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/evaluation/sim_PRS_real_MIXPRS_compare_r2.csv"))
long_table <- melt(prs_table, id.vars = c("n", "pop", "p", "rhog", "sample1", "sample2"),
                   variable.name = "method", value.name = "r2")
long_table$GWAS_type = ifelse(long_table$method == "MIXPRS_full_FALSE", "full", "prune")
long_table$approx = ifelse(long_table$method == "MIXPRS_prune_TRUE", "Identity", "Reference LD")
long_table$pop <- factor(long_table$pop, levels = c("EUR","EAS","AFR", "SAS", "AMR"))
long_table$p <- factor(long_table$p,
                       levels = c(0.1, 0.01, 0.001, 5e-04),
                       labels = c("p = 0.1", "p = 0.01", "p = 0.001", "p = 5 × 10⁻⁴"))

long_table$method <- dplyr::recode(
  long_table$method,
  "MIXPRS_full_FALSE"  = "RefLD_Full",
  "MIXPRS_prune_FALSE" = "RefLD_Prune",
  "MIXPRS_prune_TRUE"  = "Identity_Prune"
)

summary_table <- long_table %>%
  group_by(pop, p, method) %>%
  summarise(mean_r2 = mean(r2), se_r2 = sd(r2)/sqrt(n()), .groups = "drop")

p2 <- ggplot(summary_table, aes(x = pop, y = mean_r2, color = method, group = method)) +
  geom_point(size = 1, shape = 16) +
  facet_wrap(~ p, ncol = 4) +
  labs(x = NULL, y = expression(R^2), color = "Method") +
  theme_classic() +
  my_theme +
  scale_y_continuous(limits = c(NA, 0.41))

combined_plot <- ggarrange(p1, p2, ncol = 1, nrow = 2, labels = c("a","b"), heights = c(0.5,0.5))
print(combined_plot)

# Save the figure to match the journal guidelines
ggsave(filename = "FigureS2.pdf",
       plot = combined_plot,
       width = 6.5, 
       height = 8,
       units = "in",
       dpi = 300,
       device   = cairo_pdf)  # high-resolution for publication


## Aditional info
## average improvement
library(dplyr)
library(tidyr)

wide_table <- long_table %>%
  select(pop, p, method, r2, n) %>%       # drop se_r2
  pivot_wider(
    names_from  = method,
    values_from = r2
  )

avg_improvement <- wide_table %>%
  # compute the per-p differences
  mutate(
    diff_RefLD_Prune    = (RefLD_Prune    - RefLD_Full) / RefLD_Full * 100,
    diff_Identity_Prune = (Identity_Prune - RefLD_Full) / RefLD_Full * 100
  ) %>%
  # now average those diffs across the four values of p, grouped by pop
  group_by(pop) %>%
  summarise(
    avg_improvement_RefLD_Prune    = mean(diff_RefLD_Prune,    na.rm = TRUE),
    avg_improvement_Identity_Prune = mean(diff_Identity_Prune, na.rm = TRUE)
  )

print(avg_improvement)
