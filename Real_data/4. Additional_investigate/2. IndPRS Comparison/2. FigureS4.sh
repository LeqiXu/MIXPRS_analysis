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

# Fixed plot value
all_method_color = c("MIXPRS" = "#6A0DAD","Ridge" = "#E69F00", "Lasso" = "#D55E00", "Elastic Net"= "#0072B2","NNLS" = "#6A0DAD")

my_theme <- theme(
  # shrink axis titles
  axis.title.x = element_text(size = 8, face = "bold"),
  axis.title.y = element_text(size = 8, face = "bold"),
  # shrink tick labels
  axis.text.x  = element_text(size = 5, face = "bold"),
  axis.text.y  = element_text(size = 6, face = "bold"),
  # leave all other text at your preferred size
  plot.title   = element_text(size = 10, face = "bold", hjust = 0.5),
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

setwd("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/plot")

# Continuous trait
UKBB_continuous_table = fread("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/evaluation/IndPRS_Lasso_Ridge_Enet_NNLS_vs_MIXPRS_Evaluation.csv")

UKBB_continuous_table = UKBB_continuous_table[,c("population", "trait", "mixprs_r2_mean", "ridge_r2_mean","lasso_r2_mean","enet_r2_mean","nnls_r2_mean")]
colnames(UKBB_continuous_table) = c("pop", "trait", "MIXPRS", "Ridge","Lasso","Elastic Net","NNLS")

## Reshape the data to a long format and estimate the mean and sd
UKBB_long_table <- melt(UKBB_continuous_table, id.vars = c("pop", "trait"),
                   variable.name = "method", value.name = "mean_metric")
UKBB_long_table$trait = factor(UKBB_long_table$trait, levels = c("HDL","LDL","TC","logTG","Height","BMI","SBP","DBP","PLT","WBC","NEU","LYM","MON","EOS","RBC","HB","HCT","MCH","MCV","ALT","ALP","GGT","T2D","BrC","CAD","LuC"))
UKBB_long_table$pop <- factor(UKBB_long_table$pop, levels = c("EAS","AFR","SAS","AMR"))
UKBB_long_table$method = factor(UKBB_long_table$method, levels = c("MIXPRS","Ridge","Lasso","Elastic Net","NNLS"))
UKBB_long_table$mean_metric[which(UKBB_long_table$mean_metric == 0)] = NA

UKBB_long_table <- UKBB_long_table %>%
  group_by(trait, pop) %>%
  arrange(desc(mean_metric)) %>%
  mutate(rank = row_number()) %>%
  ungroup() %>%
  mutate(annotation = case_when(
    rank == 1 ~ "**",  # Best method
    rank == 2 ~ "*",   # Second best method
    TRUE ~ ""          # Others
  ))

## Obtain the legend and plot
## Calcualte the relative change of other methods over MIXPRS
UKBB_MIXPRS_ref <- UKBB_long_table[UKBB_long_table$method == "MIXPRS", c("pop", "trait", "mean_metric")]
colnames(UKBB_MIXPRS_ref)[3] <- "MIXPRS_mean_metric"
UKBB_long_table_with_MIXPRS  <- merge(UKBB_long_table, UKBB_MIXPRS_ref, by = c("pop", "trait"))

UKBB_long_table_with_MIXPRS$relative_change <- with(UKBB_long_table_with_MIXPRS , 
                                                 (mean_metric - MIXPRS_mean_metric) / MIXPRS_mean_metric)
UKBB_long_table_with_MIXPRS <- UKBB_long_table_with_MIXPRS [UKBB_long_table_with_MIXPRS $method != "MIXPRS",]
UKBB_long_table_with_MIXPRS <- UKBB_long_table_with_MIXPRS [,c("pop", "trait","method","relative_change")]
UKBB_long_table_with_MIXPRS$method = factor(UKBB_long_table_with_MIXPRS$method, levels = c("Ridge","Lasso","Elastic Net","NNLS"))

## Plot for each cohort and pop
UKBB_EAS_p = ggplot(UKBB_long_table_with_MIXPRS[which(UKBB_long_table_with_MIXPRS$pop == "EAS"),],
                  aes(x = method, y = relative_change, group = method, fill = method)) +
  geom_violin(trim=TRUE, position=position_dodge(width=0.8), width = 0.8, color="white", adjust=1.5, scale="width") +
  geom_point(aes(group=pop), position=position_jitterdodge(jitter.width = 0, dodge.width=0.5), 
             color="black", size=0.3, alpha=0.6) +
  stat_summary(fun=mean, geom="crossbar", width=0.5, 
               position=position_dodge(width=0.4), color = alpha("black", alpha=0.6), size = 0.2) +
  geom_hline(yintercept = 0, linetype = "solid", color = alpha("darkred", alpha=0.5), size = 0.6) +
  scale_fill_manual(values = all_method_color, name = "All method") +
  facet_grid(~ pop) +
  labs(title = "UKBB EAS traits",
       x = "",
       y = expression(paste("R"^2, "/", "R"^2, scriptstyle("MIXPRS"), "-1"))) +
  theme_classic() +
  my_theme +
  theme(legend.position = "none") +
  scale_y_continuous(labels = scales::percent, limits = c(-0.4,0.4))

UKBB_AFR_p = ggplot(UKBB_long_table_with_MIXPRS[which(UKBB_long_table_with_MIXPRS$pop == "AFR"),],
                  aes(x = method, y = relative_change, group = method, fill = method)) +
  geom_violin(trim=TRUE, position=position_dodge(width=0.8), width = 0.8, color="white", adjust=1.5, scale="width") +
  geom_point(aes(group=pop), position=position_jitterdodge(jitter.width = 0, dodge.width=0.5), 
             color="black", size=0.3, alpha=0.6) +
  stat_summary(fun=mean, geom="crossbar", width=0.5, 
               position=position_dodge(width=0.4), color = alpha("black", alpha=0.6), size = 0.2) +
  geom_hline(yintercept = 0, linetype = "solid", color = alpha("darkred", alpha=0.5), size = 0.6) +
  scale_fill_manual(values = all_method_color, name = "All method") +
  facet_grid(~ pop) +
  labs(title = "UKBB AFR traits",
       x = "",
       y = expression(paste("R"^2, "/", "R"^2, scriptstyle("MIXPRS"), "-1"))) +
  theme_classic() +
  my_theme +
  theme(legend.position = "none") +
  scale_y_continuous(labels = scales::percent, limits = c(-0.4,0.4))

UKBB_SAS_p = ggplot(UKBB_long_table_with_MIXPRS[which(UKBB_long_table_with_MIXPRS$pop == "SAS"),],
                  aes(x = method, y = relative_change, group = method, fill = method)) +
  geom_violin(trim=TRUE, position=position_dodge(width=0.8), width = 0.8, color="white", adjust=1.5, scale="width") +
  geom_point(aes(group=pop), position=position_jitterdodge(jitter.width = 0, dodge.width=0.5), 
             color="black", size=0.3, alpha=0.6) +
  stat_summary(fun=mean, geom="crossbar", width=0.5, 
               position=position_dodge(width=0.4), color = alpha("black", alpha=0.6), size = 0.2) +
  geom_hline(yintercept = 0, linetype = "solid", color = alpha("darkred", alpha=0.5), size = 0.6) +
  scale_fill_manual(values = all_method_color, name = "All method") +
  facet_grid(~ pop) +
  labs(title = "UKBB SAS traits",
       x = "",
       y = expression(paste("R"^2, "/", "R"^2, scriptstyle("MIXPRS"), "-1"))) +
  theme_classic() +
  my_theme +
  theme(legend.position = "none") +
  scale_y_continuous(labels = scales::percent, limits = c(-0.4,0.4))

  
UKBB_AMR_p = ggplot(UKBB_long_table_with_MIXPRS[which(UKBB_long_table_with_MIXPRS$pop == "AMR"),],
                  aes(x = method, y = relative_change, group = method, fill = method)) +
  geom_violin(trim=TRUE, position=position_dodge(width=0.8), width = 0.8, color="white", adjust=1.5, scale="width") +
  geom_point(aes(group=pop), position=position_jitterdodge(jitter.width = 0, dodge.width=0.5), 
             color="black", size=0.3, alpha=0.6) +
  stat_summary(fun=mean, geom="crossbar", width=0.5, 
               position=position_dodge(width=0.4), color = alpha("black", alpha=0.6), size = 0.2) +
  geom_hline(yintercept = 0, linetype = "solid", color = alpha("darkred", alpha=0.5), size = 0.6) +
  scale_fill_manual(values = all_method_color, name = "All method") +
  facet_grid(~ pop) +
  labs(title = "UKBB AMR traits",
       x = "",
       y = expression(paste("R"^2, "/", "R"^2, scriptstyle("MIXPRS"), "-1"))) +
  theme_classic() +
  my_theme +
  theme(legend.position = "none") +
  scale_y_continuous(labels = scales::percent, limits = c(-0.4,0.4))

combined_plot <- ggarrange(UKBB_EAS_p, UKBB_AFR_p, UKBB_SAS_p, UKBB_AMR_p, labels = c("a","b","c","d"), ncol = 2, nrow = 2)
print(combined_plot)

# Save the figure to match the journal guidelines
ggsave(filename = "FigureS4.pdf",
       plot = combined_plot,
       width = 6.5, 
       height = 8,
       units = "in",
       dpi = 300,
       device   = cairo_pdf)  # high-resolution for publication