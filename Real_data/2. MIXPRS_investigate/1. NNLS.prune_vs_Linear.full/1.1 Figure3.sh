# Step1: Data prepare
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
## prune NNLS table
r2_table_GLGC = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/evaluation/GLGC_MIX_prune_PRS_r2_Linear_vs_NNLS.csv"))
r2_table_PAGE = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/evaluation/PAGE_MIX_prune_PRS_r2_Linear_vs_NNLS.csv"))
r2_table_BBJ = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/evaluation/BBJ_MIX_prune_PRS_r2_Linear_vs_NNLS.csv"))

r2_table_GLGC = r2_table_GLGC[,c("GWAS_type","approx","trait","target_pop","snplist","selection_criterion","weight_r2")]
r2_table_PAGE = r2_table_PAGE[,c("GWAS_type","approx","trait","target_pop","snplist","selection_criterion","weight_r2")]
r2_table_BBJ = r2_table_BBJ[,c("GWAS_type","approx","trait","target_pop","snplist","selection_criterion","weight_r2")]

r2_table = rbind(r2_table_GLGC,r2_table_PAGE)
r2_table = rbind(r2_table,r2_table_BBJ)
r2_table$target_pop = factor(r2_table$target_pop, levels = c("EUR","EAS","AFR","SAS","AMR"))

r2_table_NNLS = r2_table[which(r2_table$selection_criterion == "NNLS"),c("trait","target_pop","weight_r2")]
colnames(r2_table_NNLS) = c("trait","target_pop","prune_NNLS_r2")

## full linear table
r2_table_GLGC = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/evaluation/GLGC_MIX_prune_PRS_r2_full_vs_prune.csv"))
r2_table_PAGE = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/evaluation/PAGE_MIX_prune_PRS_r2_full_vs_prune.csv"))
r2_table_BBJ = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/evaluation/BBJ_MIX_prune_PRS_r2_full_vs_prune.csv"))

r2_table_GLGC = r2_table_GLGC[,c("GWAS_type","approx","trait","target_pop","snplist","weight_r2")]
r2_table_PAGE = r2_table_PAGE[,c("GWAS_type","approx","trait","target_pop","snplist","weight_r2")]
r2_table_BBJ = r2_table_BBJ[,c("GWAS_type","approx","trait","target_pop","snplist","weight_r2")]

r2_table = rbind(r2_table_GLGC,r2_table_PAGE)
r2_table = rbind(r2_table,r2_table_BBJ)
r2_table$target_pop = factor(r2_table$target_pop, levels = c("EUR","EAS","AFR","SAS","AMR"))

r2_table_full = r2_table[which(r2_table$GWAS_type == "subsample_full"),c("trait","target_pop","weight_r2")]
colnames(r2_table_full) = c("trait","target_pop","full_Linear_r2")

## prune NNLS vs full Linear
r2_prune_NNLS_vs_full_Linear = merge(r2_table_NNLS,r2_table_full, by = c("trait","target_pop"))
r2_prune_NNLS_vs_full_Linear$prune_NNLS_improve_over_full_Linear = (r2_prune_NNLS_vs_full_Linear$prune_NNLS_r2 - r2_prune_NNLS_vs_full_Linear$full_Linear_r2) / r2_prune_NNLS_vs_full_Linear$full_Linear_r2 * 100
r2_prune_NNLS_vs_full_Linear = r2_prune_NNLS_vs_full_Linear[,c("trait","target_pop","prune_NNLS_improve_over_full_Linear")]


## AUC_table
## prune NNLS table
AUC_table_Binary_3 = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/evaluation/Binary_3_MIX_prune_PRS_AUC_Linear_vs_NNLS.csv"))
AUC_table_Binary_2 = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/evaluation/Binary_2_MIX_prune_PRS_AUC_Linear_vs_NNLS.csv"))

AUC_table_Binary_3 = AUC_table_Binary_3[,c("GWAS_type","approx","trait","target_pop","snplist","selection_criterion","weight_AUC")]
AUC_table_Binary_2 = AUC_table_Binary_2[,c("GWAS_type","approx","trait","target_pop","snplist","selection_criterion","weight_AUC")]

AUC_table = rbind(AUC_table_Binary_3,AUC_table_Binary_2)
AUC_table$target_pop = factor(AUC_table$target_pop, levels = c("EUR","EAS","AFR","SAS","AMR"))

AUC_table_NNLS = AUC_table[which(AUC_table$selection_criterion == "NNLS"),c("trait","target_pop","weight_AUC")]
colnames(AUC_table_NNLS) = c("trait","target_pop","prune_NNLS_AUC")

## full linear table
AUC_table_Binary_3 = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/evaluation/Binary_3_MIX_prune_PRS_AUC_full_vs_prune.csv"))
AUC_table_Binary_2 = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/evaluation/Binary_2_MIX_prune_PRS_AUC_full_vs_prune.csv"))

AUC_table_Binary_3 = AUC_table_Binary_3[,c("GWAS_type","approx","trait","target_pop","snplist","weight_AUC")]
AUC_table_Binary_2 = AUC_table_Binary_2[,c("GWAS_type","approx","trait","target_pop","snplist","weight_AUC")]

AUC_table = rbind(AUC_table_Binary_3,AUC_table_Binary_2)
AUC_table$target_pop = factor(AUC_table$target_pop, levels = c("EUR","EAS","AFR","SAS","AMR"))

AUC_table_full = AUC_table[which(AUC_table$GWAS_type == "subsample_full"),c("trait","target_pop","weight_AUC")]
colnames(AUC_table_full) = c("trait","target_pop","full_Linear_AUC")

## prune NNLS vs full Linear
AUC_prune_NNLS_vs_full_Linear = merge(AUC_table_NNLS,AUC_table_full, by = c("trait","target_pop"))
AUC_prune_NNLS_vs_full_Linear$prune_NNLS_improve_over_full_Linear = (AUC_prune_NNLS_vs_full_Linear$prune_NNLS_AUC - AUC_prune_NNLS_vs_full_Linear$full_Linear_AUC) / AUC_prune_NNLS_vs_full_Linear$full_Linear_AUC * 100
AUC_prune_NNLS_vs_full_Linear = AUC_prune_NNLS_vs_full_Linear[,c("trait","target_pop","prune_NNLS_improve_over_full_Linear")]


## full_table
prs_relative = rbind(r2_prune_NNLS_vs_full_Linear,AUC_prune_NNLS_vs_full_Linear)
prs_average <- prs_relative %>%
  group_by(target_pop) %>%
  summarize(mean_improvement = mean(prune_NNLS_improve_over_full_Linear, na.rm = TRUE), .groups = "drop")
prs_average$target_pop = factor(prs_average$target_pop,levels=c("EAS","AFR","SAS","AMR"))
prs_average <- prs_average %>% arrange(target_pop)
print(prs_average)


## relative plot
EAS_prs_relative = prs_relative[prs_relative$target_pop == "EAS",]
EAS_prs_relative = EAS_prs_relative %>% arrange(prune_NNLS_improve_over_full_Linear)
EAS_prs_relative$trait = factor(EAS_prs_relative$trait, levels = unique(EAS_prs_relative$trait))
EAS_p <- ggplot(EAS_prs_relative, aes(x = trait, y = prune_NNLS_improve_over_full_Linear)) +
  geom_bar(stat = "identity", 
           position = "dodge",
           aes(fill = prune_NNLS_improve_over_full_Linear < 0)) + 
  geom_text(data = prs_average[which(prs_average$target_pop=="EAS"),], aes(x = 20/3, y = 10, 
                                    label = paste("Average increase percentage:", sprintf("%.2f%%", mean_improvement))), 
            vjust = -4, size=3, 
            fontface = "bold") + 
  scale_y_continuous(labels = scales::percent_format(scale = 1, accuracy = 1),
                     name = expression(paste("R"^2, scriptstyle("MIXPRS"), "/", "R"^2, scriptstyle("Alternative"), "-1")), 
                     limits = c(-140,140)) +
  scale_fill_manual(values = c("#FB9A99", "#1F78B4")) + 
  theme_classic(base_size = 15) +
  my_theme +
  theme(legend.position="none") +
  ggtitle("EAS traits") + xlab("")

AFR_prs_relative = prs_relative[prs_relative$target_pop == "AFR",]
AFR_prs_relative = AFR_prs_relative %>% arrange(prune_NNLS_improve_over_full_Linear)
AFR_prs_relative$trait = factor(AFR_prs_relative$trait, levels = unique(AFR_prs_relative$trait))
AFR_p <- ggplot(AFR_prs_relative, aes(x = trait, y = prune_NNLS_improve_over_full_Linear)) +
  geom_bar(stat = "identity", 
           position = "dodge",
           aes(fill = prune_NNLS_improve_over_full_Linear < 0)) + 
  geom_text(data = prs_average[which(prs_average$target_pop=="AFR"),], aes(x = 9.5/3, y = 10, 
                                    label = paste("Average increase percentage:", sprintf("%.2f%%", mean_improvement))), 
            vjust = -4, size=3, 
            fontface = "bold") + 
  scale_y_continuous(labels = scales::percent_format(scale = 1, accuracy = 1),
                     name = expression(paste("R"^2, scriptstyle("MIXPRS"), "/", "R"^2, scriptstyle("Alternative"), "-1")), 
                     limits = c(-140,140)) +
  scale_fill_manual(values = c("#FB9A99", "#1F78B4")) + 
  theme_classic(base_size = 15) +
  my_theme +
  theme(legend.position="none") +
  ggtitle("AFR traits") + xlab("")

SAS_prs_relative = prs_relative[prs_relative$target_pop == "SAS",]
SAS_prs_relative = SAS_prs_relative %>% arrange(prune_NNLS_improve_over_full_Linear)
SAS_prs_relative$trait = factor(SAS_prs_relative$trait, levels = unique(SAS_prs_relative$trait))
SAS_p <- ggplot(SAS_prs_relative, aes(x = trait, y = prune_NNLS_improve_over_full_Linear)) +
  geom_bar(stat = "identity", 
           position = "dodge",
           aes(fill = prune_NNLS_improve_over_full_Linear < 0)) + 
  geom_text(data = prs_average[which(prs_average$target_pop=="SAS"),], aes(x = 8/3, y = 10, 
                                    label = paste("Average increase percentage:", sprintf("%.2f%%", mean_improvement))), 
            vjust = -4, size=3, 
            fontface = "bold") + 
  scale_y_continuous(labels = scales::percent_format(scale = 1, accuracy = 1),
                     name = expression(paste("R"^2, scriptstyle("MIXPRS"), "/", "R"^2, scriptstyle("Alternative"), "-1")), 
                     limits = c(-140,140)) +
  scale_fill_manual(values = c("#FB9A99", "#1F78B4")) + 
  theme_classic(base_size = 15) +
  my_theme +
  theme(legend.position="none") +
  ggtitle("SAS traits") + xlab("")

AMR_prs_relative = prs_relative[prs_relative$target_pop == "AMR",]
AMR_prs_relative = AMR_prs_relative %>% arrange(prune_NNLS_improve_over_full_Linear)
AMR_prs_relative$trait = factor(AMR_prs_relative$trait, levels = unique(AMR_prs_relative$trait))
AMR_p <- ggplot(AMR_prs_relative, aes(x = trait, y = prune_NNLS_improve_over_full_Linear)) +
  geom_bar(stat = "identity", 
           position = "dodge",
           aes(fill = prune_NNLS_improve_over_full_Linear < 0)) + 
  geom_text(data = prs_average[which(prs_average$target_pop=="AMR"),], aes(x = 8/3, y = 10, 
                                    label = paste("Average increase percentage:", sprintf("%.2f%%", mean_improvement))), 
            vjust = -4, size=3, 
            fontface = "bold") + 
  scale_y_continuous(labels = scales::percent_format(scale = 1, accuracy = 1),
                     name = expression(paste("R"^2, scriptstyle("MIXPRS"), "/", "R"^2, scriptstyle("Alternative"), "-1")), 
                     limits = c(-140,140)) +
  scale_fill_manual(values = c("#FB9A99", "#1F78B4")) + 
  theme_classic(base_size = 15) +
  my_theme +
  theme(legend.position="none") +
  ggtitle("AMR traits") + xlab("")

combined_plot <- ggarrange(EAS_p, AFR_p, ggarrange(SAS_p,AMR_p,ncol=2,nrow=1, widths = c(0.5,0.5),labels = c("c","d")), labels = c("a","b"), ncol = 1, nrow = 3)
print(combined_plot)

# Save the figure to match the journal guidelines
ggsave(filename = "Figure3.pdf",
       plot = combined_plot,
       width = 6.5, 
       height = 8,
       units = "in",
       dpi = 300)  # high-resolution for publication