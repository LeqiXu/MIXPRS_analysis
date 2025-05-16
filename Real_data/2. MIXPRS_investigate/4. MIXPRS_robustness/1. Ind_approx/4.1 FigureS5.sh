# Step1: Data prepare
library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)

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

## r2_table
r2_table_GLGC = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/evaluation/GLGC_MIX_prune_PRS_r2_ind_approx.csv"))
r2_table_PAGE = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/evaluation/PAGE_MIX_prune_PRS_r2_ind_approx.csv"))
r2_table_BBJ = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/evaluation/BBJ_MIX_prune_PRS_r2_ind_approx.csv"))

r2_table_GLGC = r2_table_GLGC[,c("GWAS_type","approx","trait","target_pop","snplist","weight_r2")]
r2_table_PAGE = r2_table_PAGE[,c("GWAS_type","approx","trait","target_pop","snplist","weight_r2")]
r2_table_BBJ = r2_table_BBJ[,c("GWAS_type","approx","trait","target_pop","snplist","weight_r2")]

r2_table = rbind(r2_table_GLGC,r2_table_PAGE)
r2_table = rbind(r2_table,r2_table_BBJ)

r2_table$snplist = paste0("snplist_",r2_table$snplist)
r2_table$approx = ifelse(r2_table$approx == "TRUE","Identity","Reference LD")
r2_table$target_pop = factor(r2_table$target_pop, levels = c("EUR","EAS","AFR","SAS","AMR"))

## AUC_table
AUC_table_Binary_3 = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/evaluation/Binary_3_MIX_prune_PRS_AUC_ind_approx.csv"))
AUC_table_Binary_2 = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/evaluation/Binary_2_MIX_prune_PRS_AUC_ind_approx.csv"))

AUC_table_Binary_3 = AUC_table_Binary_3[,c("GWAS_type","approx","trait","target_pop","snplist","weight_AUC")]
AUC_table_Binary_2 = AUC_table_Binary_2[,c("GWAS_type","approx","trait","target_pop","snplist","weight_AUC")]

AUC_table = rbind(AUC_table_Binary_3,AUC_table_Binary_2)

AUC_table$snplist = paste0("snplist_",AUC_table$snplist)
AUC_table$approx = ifelse(AUC_table$approx == "TRUE","Identity","Reference LD")
AUC_table$target_pop = factor(AUC_table$target_pop, levels = c("EUR","EAS","AFR","SAS","AMR"))

# Step2: plot
r2_plot = ggplot(r2_table, aes(x = approx, y = weight_r2, group = approx, fill = approx)) +
  geom_violin(trim=TRUE, position=position_dodge(width=0.8), width = 0.8, color="white", adjust=1.5, scale="width") +
  geom_point(aes(group=target_pop), position=position_jitterdodge(jitter.width = 0, dodge.width=0.5), 
             color="black", size=0.3, alpha=0.6) +
  stat_summary(fun=mean, geom="crossbar", width=0.5, 
               position=position_dodge(width=0.4), color = alpha("black", alpha=0.6), size = 0.2) +
  facet_grid(~ target_pop) +
  labs(title = "Continuous traits",
       x = "",
       y = expression(paste("R"^2))) +
  theme_classic() +
  my_theme +
  theme(legend.position = "none") +
  scale_y_continuous(limits = c(NA, 0.25))

AUC_plot = ggplot(AUC_table, aes(x = approx, y = weight_AUC, group = approx, fill = approx)) +
  geom_violin(trim=TRUE, position=position_dodge(width=0.8), width = 0.8, color="white", adjust=1.5, scale="width") +
  geom_point(aes(group=target_pop), position=position_jitterdodge(jitter.width = 0, dodge.width=0.5), 
             color="black", size=0.3, alpha=0.6) +
  stat_summary(fun=mean, geom="crossbar", width=0.5, 
               position=position_dodge(width=0.4), color = alpha("black", alpha=0.6), size = 0.2) +
  facet_grid(~ target_pop) +
  labs(title = "Binary traits",
       x = "",
       y = expression(paste("AUC"))) +
  theme_classic() +
  my_theme +
  theme(legend.position = "none") +
  scale_y_continuous(limits = c(NA, 0.7))

combined_plot <- ggarrange(r2_plot, AUC_plot, labels = c("a","b"), ncol = 1, nrow = 2)
print(combined_plot)

# Save the figure to match the journal guidelines
ggsave(filename = "FigureS5.pdf",
       plot = combined_plot,
       width = 6.5, 
       height = 8,
       units = "in",
       dpi = 300)  # high-resolution for publication
