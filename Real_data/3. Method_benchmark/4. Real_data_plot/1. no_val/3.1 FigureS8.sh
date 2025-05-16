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
"PRS-CSx" = "#006400","MUSSEL" = "#8FBC8F","PROSPER" = "#5E92F3","BridgePRS" = "#89CFF0", "PUMAS-EN" = "#bf812d", "PUMAS-EN_paper" = "#dfc27d")

auto_method_color = c("MIXPRS" = "#6A0DAD","JointPRS-auto" = "#E65475", "XPASS" = "#FF8C00", "SDPRX" = "#fec44f","PRS-CSx-auto" = "#FDCAC7", "PUMAS-EN" = "#bf812d", "PUMAS-EN_paper" = "#dfc27d")

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
GLGC_prs_table = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/evaluation/GLGC_PUMAS_EN_PRS_r2.csv"))
colnames(GLGC_prs_table) = c("pop","trait","MIXPRS","PUMAS-EN","PUMAS-EN_paper")

## Reshape the data to a long format and estimate the mean and sd
long_table <- melt(GLGC_prs_table, id.vars = c("pop", "trait"),
                   variable.name = "method", value.name = "metric")
long_table$trait = factor(long_table$trait, levels = c("HDL","LDL","TC","logTG"))
long_table$pop <- factor(long_table$pop, levels = c("EAS","AFR","SAS","AMR"))
long_table$method = factor(long_table$method, levels = c("MIXPRS","PUMAS-EN","PUMAS-EN_paper"))
long_table$metric[which(long_table$metric == 0)] = NA

long_table <- long_table %>%
  group_by(trait, pop) %>%
  arrange(desc(metric)) %>%
  mutate(rank = row_number()) %>%
  ungroup() %>%
  mutate(annotation = case_when(
    rank == 1 ~ "**",  # Best method
    rank == 2 ~ "*",   # Second best method
    TRUE ~ ""          # Others
  ))

## Obtain the legend and plot
EAS_p <- ggplot(long_table[which(long_table$pop == "EAS"),], aes(x = trait, y = metric, group = method, fill = method)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = annotation), vjust = 0, position = position_dodge(width = 0.9), color = "red", fontface = "bold", size = 2.5) +  # Bold stars
  scale_fill_manual(values = auto_method_color) +
  facet_grid(~ pop, scales = "free_x", space = "free_x") +
  labs(title = "EAS traits",
       x = "",
       y = "R²",
       fill = "Method  ") +
  theme_classic() + 
  my_theme + scale_y_continuous(limits = c(NA, 0.25))

AFR_p <- ggplot(long_table[which(long_table$pop == "AFR"),], aes(x = trait, y = metric, group = method, fill = method)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = annotation), vjust = 0, position = position_dodge(width = 0.9), color = "red", fontface = "bold", size = 2.5) +  # Bold stars
  scale_fill_manual(values = auto_method_color) +
  facet_grid(~ pop, scales = "free_x", space = "free_x") +
  labs(title = "AFR traits",
       x = "",
       y = "R²") +
  theme_classic() + 
  my_theme + scale_y_continuous(limits = c(NA, 0.25)) +
  theme(legend.position = "none")

SAS_p <- ggplot(long_table[which(long_table$pop == "SAS"),], aes(x = trait, y = metric, group = method, fill = method)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = annotation), vjust = 0, position = position_dodge(width = 0.9), color = "red", fontface = "bold", size = 2.5) +  # Bold stars
  scale_fill_manual(values = auto_method_color) +
  facet_grid(~ pop, scales = "free_x", space = "free_x") +
  labs(title = "SAS traits",
       x = "",
       y = "R²") +
  theme_classic() + 
  my_theme + scale_y_continuous(limits = c(NA, 0.25)) +
  theme(legend.position = "none")

combined_plot <- ggarrange(EAS_p, AFR_p, SAS_p, labels = c("a","b","c"), ncol = 1, nrow = 3, common.legend = TRUE)
print(combined_plot)

# Save the figure to match the journal guidelines
ggsave(filename = "FigureS8.pdf",
       plot = combined_plot,
       width = 6.5, 
       height = 8,
       units = "in",
       dpi = 300,
       device   = cairo_pdf)  # high-resolution for publication