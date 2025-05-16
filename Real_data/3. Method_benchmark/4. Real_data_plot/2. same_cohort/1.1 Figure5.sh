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
## Data preprocessing
prs_table1 = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/evaluation/PRS_update_same_cohort_r2.csv"))
prs_table2 = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/evaluation/PRS_update_same_cohort_auc.csv"))
prs_table = rbind(prs_table1,prs_table2)
prs_table = prs_table[,c("pop", "trait", "MIX" , "JointPRS_tune_max","SDPRX_auto_2","XPASS_auto_2","PRScsx_tune_max","PROSPER_tune_max","MUSSEL_tune_max","BridgePRS_tune_2")]
colnames(prs_table) = c("pop", "trait", "MIXPRS", "JointPRS","SDPRX","XPASS","PRS-CSx","PROSPER","MUSSEL","BridgePRS")

## Reshape the data to a long format and estimate the mean and sd
long_table <- melt(prs_table, id.vars = c("pop", "trait"),
                   variable.name = "method", value.name = "metric")
long_table$trait = factor(long_table$trait, levels = c("HDL","LDL","TC","logTG","Height","BMI","SBP","DBP","PLT","WBC","NEU","LYM","MON","EOS","RBC","HB","HCT","MCH","MCV","ALT","ALP","GGT","T2D","BrC","CAD","LuC"))
long_table$pop <- factor(long_table$pop, levels = c("EAS","AFR","SAS","AMR"))
long_table$method = factor(long_table$method, levels = c("MIXPRS","JointPRS","XPASS","SDPRX","PRS-CSx","MUSSEL","PROSPER","BridgePRS"))
long_table$metric[which(long_table$metric == 0)] = NA

long_table <- long_table %>%
  group_by(pop, trait, method) %>%
  summarise(
    mean_metric = mean(metric, na.rm = TRUE),
    sd_metric = sd(metric, na.rm = TRUE),
    .groups = 'drop' # This argument drops the grouping structure afterwards
  )

GLGC_trait = c("HDL","LDL","TC","logTG")
PAGE_trait = c("Height","BMI","SBP","DBP","PLT")
BBJ_trait = c("WBC","NEU","LYM","MON","EOS","RBC","HB","HCT","MCH","MCV","ALT","ALP","GGT")

## Calcualte the relative change of other methods over MIXPRS
MIXPRS_ref <- long_table[long_table$method == "MIXPRS", c("pop", "trait", "mean_metric")]
colnames(MIXPRS_ref)[3] <- "MIXPRS_mean_metric"
long_table_with_MIXPRS  <- merge(long_table, MIXPRS_ref, by = c("pop", "trait"))

long_table_with_MIXPRS$relative_change <- with(long_table_with_MIXPRS , 
                                                 (mean_metric - MIXPRS_mean_metric) / MIXPRS_mean_metric)
long_table_with_MIXPRS <- long_table_with_MIXPRS [long_table_with_MIXPRS $method != "MIXPRS",]
long_table_with_MIXPRS <- long_table_with_MIXPRS [,c("pop", "trait","method","relative_change")]
long_table_with_MIXPRS$method = factor(long_table_with_MIXPRS$method, levels = c("JointPRS", "XPASS","SDPRX","PRS-CSx","MUSSEL","PROSPER","BridgePRS"))

## Plot for each cohort and pop
EAS_p = ggplot(long_table_with_MIXPRS[which(long_table_with_MIXPRS$pop == "EAS"),],
                  aes(x = method, y = relative_change, group = method, fill = method)) +
  geom_violin(trim=TRUE, position=position_dodge(width=0.8), width = 0.8, color="white", adjust=1.5, scale="width") +
  geom_point(aes(group=pop), position=position_jitterdodge(jitter.width = 0, dodge.width=0.5), 
             color="black", size=0.3, alpha=0.6) +
  stat_summary(fun=mean, geom="crossbar", width=0.5, 
               position=position_dodge(width=0.4), color = alpha("black", alpha=0.6), size = 0.2) +
  geom_hline(yintercept = 0, linetype = "solid", color = alpha("darkred", alpha=0.5), size = 0.6) +
  scale_fill_manual(values = all_method_color, name = "All method") +
  facet_grid(~ pop) +
  labs(title = "EAS traits",
       x = "",
       y = expression(paste("R"^2, "/", "R"^2, scriptstyle("MIXPRS"), "-1"))) +
  theme_classic() +
  my_theme +
  theme(legend.position = "none") +
  scale_y_continuous(labels = scales::percent, limits = c(-1,1))

AFR_p = ggplot(long_table_with_MIXPRS[which(long_table_with_MIXPRS$pop == "AFR"),],
                  aes(x = method, y = relative_change, group = method, fill = method)) +
  geom_violin(trim=TRUE, position=position_dodge(width=0.8), width = 0.8, color="white", adjust=1.5, scale="width") +
  geom_point(aes(group=pop), position=position_jitterdodge(jitter.width = 0, dodge.width=0.5), 
             color="black", size=0.3, alpha=0.6) +
  stat_summary(fun=mean, geom="crossbar", width=0.5, 
               position=position_dodge(width=0.4), color = alpha("black", alpha=0.6), size = 0.2) +
  geom_hline(yintercept = 0, linetype = "solid", color = alpha("darkred", alpha=0.5), size = 0.6) +
  scale_fill_manual(values = all_method_color, name = "All method") +
  facet_grid(~ pop) +
  labs(title = "AFR traits",
       x = "",
       y = expression(paste("R"^2, "/", "R"^2, scriptstyle("MIXPRS"), "-1"))) +
  theme_classic() +
  my_theme +
  theme(legend.position = "none") +
  scale_y_continuous(labels = scales::percent, limits = c(-1,1))

SAS_p = ggplot(long_table_with_MIXPRS[which(long_table_with_MIXPRS$pop == "SAS"),],
                  aes(x = method, y = relative_change, group = method, fill = method)) +
  geom_violin(trim=TRUE, position=position_dodge(width=0.8), width = 0.8, color="white", adjust=1.5, scale="width") +
  geom_point(aes(group=pop), position=position_jitterdodge(jitter.width = 0, dodge.width=0.5), 
             color="black", size=0.3, alpha=0.6) +
  stat_summary(fun=mean, geom="crossbar", width=0.5, 
               position=position_dodge(width=0.4), color = alpha("black", alpha=0.6), size = 0.2) +
  geom_hline(yintercept = 0, linetype = "solid", color = alpha("darkred", alpha=0.5), size = 0.6) +
  scale_fill_manual(values = all_method_color, name = "All method") +
  facet_grid(~ pop) +
  labs(title = "SAS traits",
       x = "",
       y = expression(paste("R"^2, "/", "R"^2, scriptstyle("MIXPRS"), "-1"))) +
  theme_classic() +
  my_theme +
  theme(legend.position = "none") +
  scale_y_continuous(labels = scales::percent, limits = c(-1,1))

  
AMR_p = ggplot(long_table_with_MIXPRS[which(long_table_with_MIXPRS$pop == "AMR"),],
                  aes(x = method, y = relative_change, group = method, fill = method)) +
  geom_violin(trim=TRUE, position=position_dodge(width=0.8), width = 0.8, color="white", adjust=1.5, scale="width") +
  geom_point(aes(group=pop), position=position_jitterdodge(jitter.width = 0, dodge.width=0.5), 
             color="black", size=0.3, alpha=0.6) +
  stat_summary(fun=mean, geom="crossbar", width=0.5, 
               position=position_dodge(width=0.4), color = alpha("black", alpha=0.6), size = 0.2) +
  geom_hline(yintercept = 0, linetype = "solid", color = alpha("darkred", alpha=0.5), size = 0.6) +
  scale_fill_manual(values = all_method_color, name = "All method") +
  facet_grid(~ pop) +
  labs(title = "AMR traits",
       x = "",
       y = expression(paste("R"^2, "/", "R"^2, scriptstyle("MIXPRS"), "-1"))) +
  theme_classic() +
  my_theme +
  theme(legend.position = "none") +
  scale_y_continuous(labels = scales::percent, limits = c(-1,1))

combined_plot <- ggarrange(EAS_p, AFR_p, SAS_p, AMR_p, labels = c("a","b","c","d"), ncol = 2, nrow = 2)
print(combined_plot)


# Save the figure to match the journal guidelines
ggsave(filename = "Figure5.pdf",
       plot = combined_plot,
       width = 6.5, 
       height = 8,
       units = "in",
       dpi = 300)  # high-resolution for publication