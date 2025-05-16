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
all_method_color = c("MIXPRS" = "#6A0DAD","JointPRS" = "#B0003C",
"JointPRS-auto" = "#E65475", "XPASS" = "#FF8C00", "SDPRX" = "#fec44f","PRS-CSx-auto" = "#FDCAC7", 
"PRS-CSx" = "#006400","MUSSEL" = "#8FBC8F","PROSPER" = "#5E92F3","BridgePRS" = "#89CFF0")

auto_method_color = c("MIXPRS" = "#6A0DAD","JointPRS-auto" = "#E65475", "XPASS" = "#FF8C00", "SDPRX" = "#fec44f","PRS-CSx-auto" = "#FDCAC7")

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

setwd("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/plot")

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

## Obtain the legend and plot
EAS_df <- long_table %>%
  filter(pop == "EAS") %>%
  mutate(is_mixprs = if_else(method == "MIXPRS", 1, 0)) %>%
  arrange(is_mixprs)

EAS_p_continuous <- ggplot(EAS_df[which(EAS_df$trait_type == "Continuous"),], aes(x = trait, y = metric, color = method)) +
  geom_point(size = 1, shape = 16) +
  scale_color_manual(values = auto_method_color) +
  facet_grid(~ pop, scales = "free_x", space = "free_x") +
  labs(title = "Continuous traits", x = "", y = "R²", color = "Method") +
  theme_classic() +
  my_theme + 
  scale_y_continuous(limits = c(NA, 0.25))

AFR_df <- long_table %>%
  filter(pop == "AFR") %>%
  mutate(is_mixprs = if_else(method == "MIXPRS", 1, 0)) %>%
  arrange(is_mixprs)

AFR_p_continuous <- ggplot(AFR_df[which(AFR_df$trait_type == "Continuous"),], aes(x = trait, y = metric, color = method)) +
  geom_point(size = 1, shape = 16) +
  scale_color_manual(values = auto_method_color) +
  facet_grid(~ pop, scales = "free_x", space = "free_x") +
  labs(title = "Continuous traits", x = "", y = "R²", color = "Method") +
  theme_classic() +
  my_theme + 
  scale_y_continuous(limits = c(NA, 0.25)) +
  theme(legend.position = "none")

SAS_df <- long_table %>%
  filter(pop == "SAS") %>%
  mutate(is_mixprs = if_else(method == "MIXPRS", 1, 0)) %>%
  arrange(is_mixprs)

SAS_p <- ggplot(SAS_df, aes(x = trait, y = metric, color = method)) +
  geom_point(size = 1, shape = 16) +
  scale_color_manual(values = auto_method_color) +
  facet_grid(~ pop, scales = "free_x", space = "free_x") +
  labs(title = "Continuous traits", x = "", y = "R²", color = "Method") +
  theme_classic() +
  my_theme + 
  scale_y_continuous(limits = c(NA, 0.25))+
  theme(legend.position = "none")

AMR_df <- long_table %>%
  filter(pop == "AMR") %>%
  mutate(is_mixprs = if_else(method == "MIXPRS", 1, 0)) %>%
  arrange(is_mixprs)

AMR_p <- ggplot(AMR_df, aes(x = trait, y = metric, color = method)) +
  geom_point(size = 1, shape = 16) +
  scale_color_manual(values = auto_method_color) +
  facet_grid(~ pop, scales = "free_x", space = "free_x") +
  labs(title = "Continuous traits", x = "", y = "R²", color = "Method") +
  theme_classic() +
  my_theme + 
  scale_y_continuous(limits = c(NA, 0.25))+
  theme(legend.position = "none")

continuous_plot <- ggarrange(EAS_p_continuous, AFR_p_continuous, ggarrange(SAS_p, AMR_p, ncol=2, nrow=1, widths = c(0.5,0.5), labels = c("c","d")), labels = c("a","b"), ncol = 1, nrow = 3, common.legend = TRUE)

EAS_p_binary <- ggplot(EAS_df[which(EAS_df$trait_type == "Binary"),], aes(x = trait, y = metric, color = method)) +
  geom_point(size = 1, shape = 16) +
  scale_color_manual(values = auto_method_color) +
  facet_grid(~ pop, scales = "free_x", space = "free_x") +
  labs(title = "Binary traits", x = "", y = "AUC", color = "Method") +
  theme_classic() +
  my_theme + 
  scale_y_continuous(limits = c(NA, 0.7))+
  theme(legend.position = "none")

AFR_p_binary <- ggplot(AFR_df[which(AFR_df$trait_type == "Binary"),], aes(x = trait, y = metric, color = method)) +
  geom_point(size = 1, shape = 16) +
  scale_color_manual(values = auto_method_color) +
  facet_grid(~ pop, scales = "free_x", space = "free_x") +
  labs(title = "Binary traits", x = "", y = "AUC", color = "Method") +
  theme_classic() +
  my_theme + 
  scale_y_continuous(limits = c(NA, 0.7)) +
  theme(legend.position = "none")

binary_plot <- ggarrange(EAS_p_binary, AFR_p_binary, labels = c("e","f"), widths = c(0.5,0.5), ncol = 2, nrow = 1)

combined_plot <- ggarrange(continuous_plot,binary_plot, heights = c(0.8,0.2), ncol = 1, nrow = 2)
print(combined_plot)

# Save the figure to match the journal guidelines
ggsave(filename = "FigureS7.pdf",
       plot = combined_plot,
       width = 6.5, 
       height = 8,
       units = "in",
       dpi = 300,
       device   = cairo_pdf)  # high-resolution for publication
