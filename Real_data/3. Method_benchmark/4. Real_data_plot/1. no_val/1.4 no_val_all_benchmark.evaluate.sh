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

# Continuous trait
GLGC_prs_table = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/evaluation/GLGC_PRS_r2.csv"))
PAGE_prs_table = fread("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/evaluation/PAGE_PRS_r2.csv")
BBJ_prs_table = fread("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/evaluation/BBJ_PRS_r2.csv")

continuous_table = rbind(GLGC_prs_table,PAGE_prs_table)
continuous_table = rbind(continuous_table,BBJ_prs_table)

continuous_table = continuous_table[,c("pop", "trait", "MIX", "JointPRS_auto_max","PRScsx_auto_max","SDPRX_auto_2","XPASS_auto_2")]
colnames(continuous_table) = c("pop", "trait", "MIXPRS", "JointPRS-auto","PRS-CSx-auto","SDPRX","XPASS")
continuous_table$trait_type = "Continuous"

# Binary trait
Binary_3_prs_table = fread("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/evaluation/Binary_3_PRS_AUC.csv")
Binary_2_prs_table = fread("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/evaluation/Binary_2_PRS_AUC.csv")

binary_table = rbind(Binary_3_prs_table,Binary_2_prs_table)

binary_table = binary_table[,c("pop", "trait", "MIX", "JointPRS_auto_max","PRScsx_auto_max","SDPRX_auto_2","XPASS_auto_2")]
colnames(binary_table) = c("pop", "trait", "MIXPRS", "JointPRS-auto","PRS-CSx-auto","SDPRX","XPASS")
binary_table$trait_type = "Binary"

## All traits
all_table <- rbind(continuous_table,binary_table)

## Reshape the data to a long format and estimate the mean and sd
long_table <- melt(all_table, id.vars = c("pop", "trait", "trait_type"),
                   variable.name = "method", value.name = "metric")
long_table$trait_type = factor(long_table$trait_type, levels = c("Continuous","Binary"))
long_table$trait = factor(long_table$trait, levels = c("HDL","LDL","TC","logTG","Height","BMI","SBP","DBP","PLT","WBC","NEU","LYM","MON","EOS","RBC","HB","HCT","MCH","MCV","ALT","ALP","GGT","T2D","BrC","CAD","LuC"))
long_table$pop <- factor(long_table$pop, levels = c("EAS","AFR","SAS","AMR"))
long_table$method = factor(long_table$method, levels = c("MIXPRS","JointPRS-auto","SDPRX","PRS-CSx-auto","XPASS"))
long_table$metric[which(long_table$metric == 0)] = NA

# Calculate improvements relative to MIXPRS
improvement_summary <- long_table %>%
  pivot_wider(names_from = method, values_from = metric) %>%
  pivot_longer(cols = c("JointPRS-auto", "SDPRX", "PRS-CSx-auto", "XPASS"),
               names_to = "Compared_Method",
               values_to = "Compared_Metric") %>%
  mutate(Improvement = (MIXPRS - Compared_Metric) / Compared_Metric) %>%
  group_by(pop, Compared_Method) %>%
  summarize(
    avg_improvement = mean(Improvement, na.rm = TRUE) * 100,
    count = sum(!is.na(Improvement)),
    .groups = "drop"
  )

improvement_summary$avg_improvement <- round(improvement_summary$avg_improvement,2)
# Display the summarized improvement table
print(improvement_summary)
