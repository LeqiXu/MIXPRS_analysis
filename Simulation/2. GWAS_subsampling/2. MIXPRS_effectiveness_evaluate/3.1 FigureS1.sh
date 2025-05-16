## MIXPRS comparison
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
library(ggsci)

setwd("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/plot")

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

auto_method_color = c("MIXPRS" = "#6A0DAD","JointPRS-auto" = "#E65475", "XPASS" = "#FF8C00", "SDPRX" = "#fec44f","PRS-CSx-auto" = "#FDCAC7")

# Data preprocessing
prs_table = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/evaluation/sim_PRS_real_MIXPRS_benefit_r2.csv"))
long_table <- melt(prs_table, id.vars = c("n", "pop", "p", "rhog", "sample1", "sample2"),
                   variable.name = "method", value.name = "r2")
long_table$pop <- factor(long_table$pop, levels = c("EUR","EAS","AFR","SAS","AMR"))
long_table$p <- factor(long_table$p,
                       levels = c(0.1, 0.01, 0.001, 5e-04),
                       labels = c("p = 0.1", "p = 0.01", "p = 0.001", "p = 5 × 10⁻⁴"))

# Obtain the plot
long_table$method <- dplyr::recode(
  long_table$method,
  "JointPRS_auto"  = "JointPRS-auto",
  "SDPRX" = "SDPRX",
  "MIXPRS_prune_TRUE"  = "MIXPRS"
)

summary_table <- long_table %>%
  group_by(pop, p, method) %>%
  summarise(mean_r2 = mean(r2), se_r2 = sd(r2)/sqrt(n()), .groups = "drop")

p1 <- ggplot(summary_table, aes(x = pop, y = mean_r2, color = method, group = method)) +
  geom_point(size = 1) +
  facet_wrap(~ p, ncol = 5) +
  labs(title = "",
       x = "",
       y = "R²",
       color = "Method  ") +
  theme_classic() +
  my_theme +
  scale_color_manual(values = auto_method_color) +
  scale_y_continuous(limits = c(NA, 0.41))

summary_table <- summary_table %>%
  group_by(pop, p) %>%
  arrange(desc(mean_r2)) %>%
  mutate(rank = row_number()) %>%
  ungroup() %>%
  mutate(annotation = case_when(
    rank == 1 ~ "**",  # Best method
    rank == 2 ~ "*",   # Second best method
    TRUE ~ ""          # Others
  ))

p2 <- ggplot(summary_table, aes(x = pop, y = mean_r2, group = method, fill = method)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = annotation), vjust = 0, position = position_dodge(width = 0.9), color = "red", fontface = "bold", size = 2.5) +  # Bold stars
  scale_fill_manual(values = auto_method_color) +
  facet_grid(~ p, scales = "free_x", space = "free_x") +
  labs(title = "",
       x = "",
       y = "R²",
       fill = "Method  ") +
  theme_classic() + 
  my_theme + 
  scale_y_continuous(limits = c(NA, 0.41))

combined_plot <- ggarrange(p1, p2, ncol = 1, nrow = 2, labels = c("a","b"), heights = c(0.5,0.5))

# Save the figure to match the journal guidelines
ggsave(filename = "FigureS1.pdf",
       plot = combined_plot,
       width = 6.5, 
       height = 8,
       units = "in",
       dpi = 300,
       device   = cairo_pdf)  # high-resolution for publication


grid.text("**", 
          x = unit(0.5, "npc"), y = unit(0.82, "npc"), 
          gp = gpar(col = "red", fontsize = 14))

grid.text("Best Method", 
          x = unit(0.55, "npc"), y = unit(0.82, "npc"), 
          gp = gpar(col = "black", fontsize = 14))

grid.text("*", 
          x = unit(0.63, "npc"), y = unit(0.82, "npc"), 
          gp = gpar(col = "red", fontsize = 14))

grid.text("Second Best Method", 
          x = unit(0.7, "npc"), y = unit(0.82, "npc"), 
          gp = gpar(col = "black", fontsize = 14))
