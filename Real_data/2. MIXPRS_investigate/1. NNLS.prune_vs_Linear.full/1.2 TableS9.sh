## MIXPRS comparison
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(writexl)

setwd("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/table")

## Data preprocessing
## r2_table
## prune table
r2_table_GLGC = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/evaluation/GLGC_MIX_prune_PRS_r2_Linear_vs_NNLS.csv"))
r2_table_PAGE = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/evaluation/PAGE_MIX_prune_PRS_r2_Linear_vs_NNLS.csv"))
r2_table_BBJ = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/evaluation/BBJ_MIX_prune_PRS_r2_Linear_vs_NNLS.csv"))

r2_table_GLGC = r2_table_GLGC[,c("GWAS_type","approx","trait","target_pop","snplist","selection_criterion","weight_r2")]
r2_table_PAGE = r2_table_PAGE[,c("GWAS_type","approx","trait","target_pop","snplist","selection_criterion","weight_r2")]
r2_table_BBJ = r2_table_BBJ[,c("GWAS_type","approx","trait","target_pop","snplist","selection_criterion","weight_r2")]

r2_table = rbind(r2_table_GLGC,r2_table_PAGE)
r2_table = rbind(r2_table,r2_table_BBJ)
r2_table$target_pop = factor(r2_table$target_pop, levels = c("EUR","EAS","AFR","SAS","AMR"))
r2_table$method = paste0(r2_table$GWAS_type,"_",r2_table$selection_criterion)

r2_table_prune = r2_table[which(r2_table$selection_criterion == "NNLS" | r2_table$selection_criterion == "NO"),c("trait","method","target_pop","weight_r2")]
r2_table_prune$method = str_replace(r2_table_prune$method,"subsample_prune","Prune")
r2_table_prune$method = str_replace(r2_table_prune$method,"NO","Linear")
colnames(r2_table_prune) = c("Trait","Strategy","Pop","Metric")

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
r2_table_full$method = "Full_Linear"
colnames(r2_table_full) = c("Trait","Pop","Metric","Strategy")
r2_table_full = r2_table_full[,c("Trait","Strategy","Pop","Metric")]
r2_table = rbind(r2_table_prune,r2_table_full)

## AUC_table
## prune NNLS table
AUC_table_Binary_3 = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/evaluation/Binary_3_MIX_prune_PRS_AUC_Linear_vs_NNLS.csv"))
AUC_table_Binary_2 = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/evaluation/Binary_2_MIX_prune_PRS_AUC_Linear_vs_NNLS.csv"))

AUC_table_Binary_3 = AUC_table_Binary_3[,c("GWAS_type","approx","trait","target_pop","snplist","selection_criterion","weight_AUC")]
AUC_table_Binary_2 = AUC_table_Binary_2[,c("GWAS_type","approx","trait","target_pop","snplist","selection_criterion","weight_AUC")]

AUC_table = rbind(AUC_table_Binary_3,AUC_table_Binary_2)
AUC_table$target_pop = factor(AUC_table$target_pop, levels = c("EUR","EAS","AFR","SAS","AMR"))
AUC_table$method = paste0(AUC_table$GWAS_type,"_",AUC_table$selection_criterion)

AUC_table_prune = AUC_table[which(AUC_table$selection_criterion == "NNLS" | AUC_table$selection_criterion == "NO"),c("trait","method","target_pop","weight_AUC")]
AUC_table_prune$method = str_replace(AUC_table_prune$method,"subsample_prune","Prune")
AUC_table_prune$method = str_replace(AUC_table_prune$method,"NO","Linear")
colnames(AUC_table_prune) = c("Trait","Strategy","Pop","Metric")

## full linear table
AUC_table_Binary_3 = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/evaluation/Binary_3_MIX_prune_PRS_AUC_full_vs_prune.csv"))
AUC_table_Binary_2 = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/evaluation/Binary_2_MIX_prune_PRS_AUC_full_vs_prune.csv"))

AUC_table_Binary_3 = AUC_table_Binary_3[,c("GWAS_type","approx","trait","target_pop","snplist","weight_AUC")]
AUC_table_Binary_2 = AUC_table_Binary_2[,c("GWAS_type","approx","trait","target_pop","snplist","weight_AUC")]

AUC_table = rbind(AUC_table_Binary_3,AUC_table_Binary_2)
AUC_table$target_pop = factor(AUC_table$target_pop, levels = c("EUR","EAS","AFR","SAS","AMR"))

AUC_table_full = AUC_table[which(AUC_table$GWAS_type == "subsample_full"),c("trait","target_pop","weight_AUC")]
AUC_table_full$method = "Full_Linear"
colnames(AUC_table_full) = c("Trait","Pop","Metric","Strategy")
AUC_table_full = AUC_table_full[,c("Trait","Strategy","Pop","Metric")]
AUC_table = rbind(AUC_table_prune,AUC_table_full)

## organize the final table
long_table <- rbind(r2_table,AUC_table)
long_table$Strategy <- factor(long_table$Strategy, levels = c("Full_Linear","Prune_Linear","Prune_NNLS"))
long_table$Pop <- factor(long_table$Pop, levels = c("EUR","EAS","AFR","SAS","AMR"))
long_table$Trait <- factor(long_table$Trait,
                           levels = c("HDL","LDL","TC","logTG","Height","BMI","SBP","DBP","PLT","WBC","NEU","LYM","MON","EOS","RBC","HB","HCT","MCH","MCV","ALT","ALP","GGT","T2D","BrC","CAD","LuC"))

# Round Metric to three decimals and reshape
reshaped_table <- long_table %>%
  mutate(Metric = round(Metric, 3)) %>%
  select(Pop, Trait, Strategy, Metric) %>%
  pivot_wider(names_from = Strategy, values_from = Metric) %>%
  arrange(Pop,Trait)

# Check results
print(reshaped_table)

write_xlsx(reshaped_table, "TableS9.xlsx")