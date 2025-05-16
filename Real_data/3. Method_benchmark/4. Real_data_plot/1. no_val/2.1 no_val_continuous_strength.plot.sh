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

## Reshape the data to a long format and estimate the mean and sd
long_table <- melt(continuous_table, id.vars = c("pop", "trait"),
                   variable.name = "method", value.name = "r2")
long_table$trait = factor(long_table$trait, levels = c("HDL","LDL","TC","logTG","Height","BMI","SBP","DBP","PLT","WBC","NEU","LYM","MON","EOS","RBC","HB","HCT","MCH","MCV","ALT","ALP","GGT"))
long_table$pop <- factor(long_table$pop, levels = c("EAS","AFR","SAS","AMR"))
long_table$method = factor(long_table$method, levels = c("MIXPRS","JointPRS-auto","SDPRX"))
long_table$r2[which(long_table$r2 == 0)] = NA

## Obtain the legend and plot
EAS_df <- long_table %>%
  filter(pop == "EAS") %>%
  mutate(is_mixprs = if_else(method == "MIXPRS", 1, 0)) %>%
  arrange(is_mixprs)

EAS_p <- ggplot(EAS_df, aes(x = trait, y = r2, color = method)) +
  geom_point(position = position_dodge(width = 0), size = 3) +
  scale_color_manual(values = auto_method_color) +
  facet_grid(~ pop, scales = "free_x", space = "free_x") +
  labs(x = "", y = "R²", color = "Method") +
  theme_classic() +
  my_theme + 
  scale_y_continuous(limits = c(NA, 0.25))

AFR_df <- long_table %>%
  filter(pop == "AFR") %>%
  mutate(is_mixprs = if_else(method == "MIXPRS", 1, 0)) %>%
  arrange(is_mixprs)

AFR_p <- ggplot(AFR_df, aes(x = trait, y = r2, color = method)) +
  geom_point(position = position_dodge(width = 0), size = 3) +
  scale_color_manual(values = auto_method_color) +
  facet_grid(~ pop, scales = "free_x", space = "free_x") +
  labs(x = "", y = "R²", color = "Method") +
  theme_classic() +
  my_theme + 
  scale_y_continuous(limits = c(NA, 0.25)) +
  no_legend

SAS_df <- long_table %>%
  filter(pop == "SAS") %>%
  mutate(is_mixprs = if_else(method == "MIXPRS", 1, 0)) %>%
  arrange(is_mixprs)

SAS_p <- ggplot(SAS_df, aes(x = trait, y = r2, color = method)) +
  geom_point(position = position_dodge(width = 0), size = 3) +
  scale_color_manual(values = auto_method_color) +
  facet_grid(~ pop, scales = "free_x", space = "free_x") +
  labs(x = "", y = "R²", color = "Method") +
  theme_classic() +
  my_theme + 
  scale_y_continuous(limits = c(NA, 0.25))+
  no_legend

AMR_df <- long_table %>%
  filter(pop == "AMR") %>%
  mutate(is_mixprs = if_else(method == "MIXPRS", 1, 0)) %>%
  arrange(is_mixprs)

AMR_p <- ggplot(AMR_df, aes(x = trait, y = r2, color = method)) +
  geom_point(position = position_dodge(width = 0), size = 3) +
  scale_color_manual(values = auto_method_color) +
  facet_grid(~ pop, scales = "free_x", space = "free_x") +
  labs(x = "", y = "R²", color = "Method") +
  theme_classic() +
  my_theme + 
  scale_y_continuous(limits = c(NA, 0.25))+
  no_legend

combined_plot <- ggarrange(EAS_p, AFR_p, ggarrange(SAS_p, AMR_p, ncol=2, nrow=1, widths = c(0.5,0.5), labels = c("c","d")), labels = c("a","b"), ncol = 1, nrow = 3, common.legend = TRUE)
print(combined_plot)

# no_val_continuous_strength.png width 1600 height 1600
