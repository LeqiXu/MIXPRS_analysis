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

setwd("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/plot")

# Fixed plot value
all_method_color = c("MIXPRS" = "#6A0DAD","JointPRS" = "#B0003C",
"JointPRS-auto" = "#E65475", "XPASS" = "#FF8C00", "SDPRX" = "#fec44f","PRS-CSx-auto" = "#FDCAC7", 
"PRS-CSx" = "#006400","MUSSEL" = "#8FBC8F","PROSPER" = "#5E92F3","BridgePRS" = "#89CFF0")

auto_method_color = c("MIXPRS" = "#6A0DAD","JointPRS-auto" = "#E65475", "XPASS" = "#FF8C00", "SDPRX" = "#fec44f","PRS-CSx-auto" = "#FDCAC7")

my_theme <- theme(
    plot.title = element_text(size=16, face = "bold"),
    text = element_text(size=16),
    axis.text.x = element_text(size=14),
    axis.text.y = element_text(size=14),
    axis.title.x = element_text(size=16, face = "bold"),
    axis.title.y = element_text(size=16),
    strip.text = element_text(face = "bold", size = 14),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey", size = 0.5),
    panel.grid.minor.y = element_line(color = "lightgrey", size = 0.25)
  )

have_legend <- theme(
    legend.text = element_text(size=14),
    legend.title = element_text(size=16, face = "bold"),
    legend.position = "top")

no_legend <- theme(
    legend.position = "none"
)

# Continuous trait
GLGC_prs_table = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/evaluation/GLGC_PRS_r2.csv"))
PAGE_prs_table = fread("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/evaluation/PAGE_PRS_r2.csv")
BBJ_prs_table = fread("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/evaluation/BBJ_PRS_r2.csv")

continuous_table = rbind(GLGC_prs_table,PAGE_prs_table)
continuous_table = rbind(continuous_table,BBJ_prs_table)

continuous_table = continuous_table[,c("pop", "trait", "MIX", "JointPRS_auto_max","SDPRX_auto_2")]
colnames(continuous_table) = c("pop", "trait", "MIXPRS", "JointPRS-auto","SDPRX")
continuous_table$trait_type = "Continuous"

# Binary trait
Binary_3_prs_table = fread("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/evaluation/Binary_3_PRS_AUC.csv")
Binary_2_prs_table = fread("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/evaluation/Binary_2_PRS_AUC.csv")

binary_table = rbind(Binary_3_prs_table,Binary_2_prs_table)

binary_table = binary_table[,c("pop", "trait", "MIX", "JointPRS_auto_max","SDPRX_auto_2")]
colnames(binary_table) = c("pop", "trait", "MIXPRS", "JointPRS-auto","SDPRX")
binary_table$trait_type = "Binary"

# All trait
all_table <- rbind(continuous_table,binary_table)

## Reshape the data to a long format and estimate the mean and sd
long_table <- melt(all_table, id.vars = c("pop", "trait", "trait_type"),
                   variable.name = "method", value.name = "metric")
long_table$trait_type = factor(long_table$trait_type, levels = c("Continuous","Binary"))
long_table$trait = factor(long_table$trait, levels = c("HDL","LDL","TC","logTG","Height","BMI","SBP","DBP","PLT","WBC","NEU","LYM","MON","EOS","RBC","HB","HCT","MCH","MCV","ALT","ALP","GGT","T2D","BrC","CAD","LuC"))
long_table$pop <- factor(long_table$pop, levels = c("EAS","AFR","SAS","AMR"))
long_table$method = factor(long_table$method, levels = c("MIXPRS","JointPRS-auto","SDPRX"))
long_table$metric[which(long_table$metric == 0)] = NA

# Count how many times each method performs best (maximum metric)
best_count <- long_table %>%
  group_by(pop, trait, trait_type) %>%
  filter(metric == max(metric, na.rm = TRUE)) %>%
  ungroup() %>%
  count(pop, trait_type, method, name = "num_best_traits") %>%
  arrange(pop, trait_type, method)

# Display the summary
print(best_count)