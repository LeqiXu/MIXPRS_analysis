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

setwd("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/plot")

# Fixed plot value
all_method_color = c("MIXPRS" = "#6A0DAD","JointPRS" = "#B0003C",
"JointPRS-auto" = "#E65475", "XPASS" = "#FF8C00", "SDPRX" = "#fec44f","PRS-CSx-auto" = "#FDCAC7", 
"PRS-CSx" = "#006400","MUSSEL" = "#8FBC8F","PROSPER" = "#5E92F3","BridgePRS" = "#89CFF0")

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

# Continuous plot
# Continuous plot
## Data preprocessing
prs_table = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/evaluation/sim_PRS_real_tune_benchmark_r2.csv"))
prs_table = prs_table[,c("n","pop","p","rhog","sample2","MIXPRS_auto_5","SDPRX_auto_2","XPASS_auto_2","JointPRS_tune_5","PRScsx_tune_5","PROSPER_tune_5","MUSSEL_tune_5","BridgePRS_tune_2")]
colnames(prs_table) = c("n","pop","p","rhog","sample2","MIXPRS","SDPRX","XPASS","JointPRS","PRS-CSx","PROSPER","MUSSEL","BridgePRS")
prs_table = prs_table %>% 
  mutate(minor_sample = case_when(sample2 == "25K" ~ "25,000",
                                  sample2 == "90K" ~ "90,000"))

## Reshape the data to a long format and estimate the mean and sd
long_table <- melt(prs_table, id.vars = c("n","pop","p","rhog","sample2","minor_sample"),
                   variable.name = "method", value.name = "r2")

long_table$p <- factor(long_table$p,
                       levels = c(0.1, 0.01, 0.001, 5e-04),
                       labels = c("p = 0.1", "p = 0.01", "p = 0.001", "p = 5 × 10⁻⁴"))
long_table$pop <- factor(long_table$pop, levels = c("EAS", "AFR", "SAS", "AMR"))
long_table$method = factor(long_table$method, levels = c("MIXPRS","JointPRS","XPASS","SDPRX","PRS-CSx","MUSSEL","PROSPER","BridgePRS"))

## Calcualte the relative change of other methods over JointPRS
long_table <- long_table %>%
  group_by(pop,p,rhog,sample2,minor_sample,method) %>%
  summarise(
    mean_r2 = mean(r2, na.rm = TRUE),
    sd_r2 = sd(r2, na.rm = TRUE),
    .groups = 'drop' # This argument drops the grouping structure afterwards
  )

MIXPRS_ref <- long_table[long_table$method == "MIXPRS", c("pop", "p", "rhog", "sample2", "minor_sample", "mean_r2")]
colnames(MIXPRS_ref)[6] <- "MIXPRS_mean_r2"
long_table_with_MIXPRS <- merge(long_table, MIXPRS_ref, by = c("pop", "p", "rhog", "sample2", "minor_sample"))

long_table_with_MIXPRS$relative_change <- with(long_table_with_MIXPRS, 
                                                 (mean_r2 - MIXPRS_mean_r2) / MIXPRS_mean_r2)
long_table_with_MIXPRS <- long_table_with_MIXPRS[long_table_with_MIXPRS$method != "MIXPRS",]
long_table_with_MIXPRS <- long_table_with_MIXPRS[,c("pop", "p", "rhog", "sample2", "minor_sample", "method", "relative_change")]
long_table_with_MIXPRS$method = factor(long_table_with_MIXPRS$method, levels = c("JointPRS","XPASS","SDPRX","PRS-CSx","MUSSEL","PROSPER","BridgePRS"))


## Plot for each cohort and pop
p_25K = ggplot(long_table_with_MIXPRS[which(long_table_with_MIXPRS$minor_sample == "25,000"),],
                  aes(x = method, y = relative_change, group = method, fill = method)) +
  geom_violin(trim=TRUE, position=position_dodge(width=0.8), width = 0.8, color="white", adjust=1.5, scale="width") +
  geom_point(aes(group=pop), position=position_jitterdodge(jitter.width = 0, dodge.width=0.5), 
             color="black", size=0.3, alpha=0.6) +
  stat_summary(fun=mean, geom="crossbar", width=0.5, 
               position=position_dodge(width=0.4), color = alpha("black", alpha=0.6), size = 0.2) +
  geom_hline(yintercept = 0, linetype = "solid", color = alpha("darkred", alpha=0.5), size = 0.6) +
  scale_fill_manual(values = all_method_color, name = "All method") +
  facet_wrap(~ pop, ncol = 2) +
  labs(title = "Non-European Population Training Sample Size = 25,000",
       x = "Method",
       y = expression(paste("R"^2, "/", "R"^2, scriptstyle("MIXPRS"), "-1"))) +
  theme_classic() +
  my_theme +
  theme(legend.position = "none") +
  scale_y_continuous(labels = scales::percent, limits = c(-1,0.2))

p_90K = ggplot(long_table_with_MIXPRS[which(long_table_with_MIXPRS$minor_sample == "90,000"),],
                  aes(x = method, y = relative_change, group = method, fill = method)) +
  geom_violin(trim=TRUE, position=position_dodge(width=0.8), width = 0.8, color="white", adjust=1.5, scale="width") +
  geom_point(aes(group=pop), position=position_jitterdodge(jitter.width = 0, dodge.width=0.5), 
             color="black", size=0.3, alpha=0.6) +
  stat_summary(fun=mean, geom="crossbar", width=0.5, 
               position=position_dodge(width=0.4), color = alpha("black", alpha=0.6), size = 0.2) +
  geom_hline(yintercept = 0, linetype = "solid", color = alpha("darkred", alpha=0.5), size = 0.6) +
  scale_fill_manual(values = all_method_color, name = "All method") +
  facet_wrap(~ pop, ncol = 2) +
  labs(title = "Non-European Population Training Sample Size = 90,000",
       x = "Method",
       y = expression(paste("R"^2, "/", "R"^2, scriptstyle("MIXPRS"), "-1"))) +
  theme_classic() +
  my_theme +
  theme(legend.position = "none") +
  scale_y_continuous(labels = scales::percent, limits = c(-1,0.2))

combined_plot <- ggarrange(p_25K, p_90K, ncol = 1, nrow = 2, labels = c("a","b"), heights = c(0.5,0.5))

# Save the figure to match the journal guidelines
ggsave(filename = "Figure2.pdf",
       plot = combined_plot,
       width = 6.5, 
       height = 8,
       units = "in",
       dpi = 300)  # high-resolution for publication
