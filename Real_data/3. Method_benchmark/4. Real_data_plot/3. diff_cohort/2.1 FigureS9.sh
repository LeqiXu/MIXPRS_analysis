## auto and tuning method comparison
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
library(ggpattern)

setwd("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/plot")

# Fixed plot value
compare_method_color = c("MIXPRS" = "#6A0DAD","IndPRS" = "#0072B2","SDPRX" = "#fec44f")

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

# Continuous trait
continuous_table = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/evaluation/MIXPRS_update_diff_cohort_r2.csv"))
continuous_table = continuous_table[,c("pop", "trait", "MIXPRS", "LinearPRS")]
colnames(continuous_table) = c("pop", "trait", "MIXPRS", "IndPRS")

prs_table = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/evaluation/PRS_update_diff_cohort_r2.csv"))
prs_table = prs_table[,c("pop", "trait", "SDPRX_auto_2")]
colnames(prs_table) = c("pop", "trait", "SDPRX")

continuous_table = merge(continuous_table,prs_table,by = c("pop", "trait"))

## Reshape the data to a long format and estimate the mean and sd
long_table <- melt(continuous_table, id.vars = c("pop", "trait"),
                   variable.name = "method", value.name = "r2")
long_table$trait = factor(long_table$trait, levels = c("HDL","LDL","TC","logTG","Height","BMI","SBP","DBP","PLT","WBC","NEU","LYM","MON","EOS","RBC","HB","HCT","MCH","MCV","ALT","ALP","GGT"))
long_table$pop <- factor(long_table$pop, levels = c("EAS","AFR","SAS","AMR"))
long_table$method = factor(long_table$method, levels = c("MIXPRS","IndPRS", "SDPRX"))
long_table$r2[which(long_table$r2 == 0)] = NA

## Obtain the legend and plot
AFR_df <- long_table %>%
  filter(pop == "AFR") %>%
  mutate(is_mixprs = if_else(method == "MIXPRS", 1, 0)) %>%
  arrange(is_mixprs)

AFR_p <- ggplot(AFR_df, aes(x = trait, y = r2, color = method)) +
  geom_point(size = 1, shape = 16) +
  scale_color_manual(values = compare_method_color) +
  facet_grid(~ pop, scales = "free_x", space = "free_x") +
  labs(title = "AFR traits", x = "", y = "R²", color = "Method") +
  theme_classic() +
  my_theme + 
  scale_y_continuous(limits = c(NA, 0.25)) +
  theme(legend.position = "none")

AMR_df <- long_table %>%
  filter(pop == "AMR") %>%
  mutate(is_mixprs = if_else(method == "MIXPRS", 1, 0)) %>%
  arrange(is_mixprs)

AMR_p <- ggplot(AMR_df, aes(x = trait, y = r2, color = method)) +
  geom_point(size = 1, shape = 16) +
  scale_color_manual(values = compare_method_color) +
  facet_grid(~ pop, scales = "free_x", space = "free_x") +
  labs(title = "AMR traits", x = "", y = "R²", color = "Method") +
  theme_classic() +
  my_theme + 
  scale_y_continuous(limits = c(NA, 0.25))+
  theme(legend.position = "none")

combined_plot <- ggarrange(AFR_p, AMR_p, labels = c("a","b"), ncol = 1, nrow = 2, common.legend = TRUE)
print(combined_plot)


# Save the figure to match the journal guidelines
ggsave(filename = "FigureS9.pdf",
       plot = combined_plot,
       width = 6.5, 
       height = 8,
       units = "in",
       dpi = 300,
       device   = cairo_pdf)  # high-resolution for publication