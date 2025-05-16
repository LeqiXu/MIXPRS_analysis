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

# Binary trait
Binary_3_prs_table = fread("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/evaluation/Binary_3_PRS_AUC.csv")
Binary_2_prs_table = fread("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/evaluation/Binary_2_PRS_AUC.csv")

binary_table = rbind(Binary_3_prs_table,Binary_2_prs_table)

binary_table = binary_table[,c("pop", "trait", "MIX", "JointPRS_auto_max","SDPRX_auto_2")]
colnames(binary_table) = c("pop", "trait", "MIXPRS", "JointPRS-auto","SDPRX")

## Reshape the data to a long format and estimate the mean and sd
long_table <- melt(binary_table, id.vars = c("pop", "trait"),
                   variable.name = "method", value.name = "AUC")
long_table$trait = factor(long_table$trait, levels = c("T2D","BrC","CAD","LuC"))
long_table$pop <- factor(long_table$pop, levels = c("EAS","AFR","SAS","AMR"))
long_table$method = factor(long_table$method, levels = c("MIXPRS","JointPRS-auto","SDPRX"))
long_table$AUC[which(long_table$AUC == 0)] = NA

## Obtain the legend and plot
EAS_df <- long_table %>%
  filter(pop == "EAS") %>%
  mutate(is_mixprs = if_else(method == "MIXPRS", 1, 0)) %>%
  arrange(is_mixprs)

EAS_p <- ggplot(EAS_df, aes(x = trait, y = AUC, color = method)) +
  geom_point(position = position_dodge(width = 0), size = 3) +
  scale_color_manual(values = auto_method_color) +
  facet_grid(~ pop, scales = "free_x", space = "free_x") +
  labs(x = "", y = "AUC", color = "Method") +
  theme_classic() +
  my_theme + 
  scale_y_continuous(limits = c(NA, 0.7))

AFR_df <- long_table %>%
  filter(pop == "AFR") %>%
  mutate(is_mixprs = if_else(method == "MIXPRS", 1, 0)) %>%
  arrange(is_mixprs)

AFR_p <- ggplot(AFR_df, aes(x = trait, y = AUC, color = method)) +
  geom_point(position = position_dodge(width = 0), size = 3) +
  scale_color_manual(values = auto_method_color) +
  facet_grid(~ pop, scales = "free_x", space = "free_x") +
  labs(x = "", y = "AUC", color = "Method") +
  theme_classic() +
  my_theme + 
  scale_y_continuous(limits = c(NA, 0.7)) +
  no_legend

combined_plot <- ggarrange(EAS_p, AFR_p, labels = c("a","b"), ncol = 1, nrow = 2, common.legend = TRUE)
print(combined_plot)

# no_val_binary_strength.png width 800 height 800
