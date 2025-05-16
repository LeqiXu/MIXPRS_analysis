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
r2_table_GLGC = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/evaluation/GLGC_MIX_prune_PRS_r2_diff_snplists.csv"))
r2_table_PAGE = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/evaluation/PAGE_MIX_prune_PRS_r2_diff_snplists.csv"))
r2_table_BBJ = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/evaluation/BBJ_MIX_prune_PRS_r2_diff_snplists.csv"))

r2_table_GLGC = r2_table_GLGC[,c("GWAS_type","approx","trait","target_pop","snplist","weight_r2")]
r2_table_PAGE = r2_table_PAGE[,c("GWAS_type","approx","trait","target_pop","snplist","weight_r2")]
r2_table_BBJ = r2_table_BBJ[,c("GWAS_type","approx","trait","target_pop","snplist","weight_r2")]

r2_table = rbind(r2_table_GLGC,r2_table_PAGE)
r2_table = rbind(r2_table,r2_table_BBJ)

r2_table$snplist = paste0("snplist_",r2_table$snplist)
r2_table$approx = ifelse(r2_table$approx == "TRUE","Identity","Reference LD")
r2_table$target_pop = factor(r2_table$target_pop, levels = c("EUR","EAS","AFR","SAS","AMR"))
r2_table = r2_table[,c("target_pop","trait","snplist","weight_r2")]
colnames(r2_table) = c("Pop","Trait","Strategy","Metric")

## AUC_table
## prune NNLS table
AUC_table_Binary_3 = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/evaluation/Binary_3_MIX_prune_PRS_AUC_diff_snplists.csv"))
AUC_table_Binary_2 = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/evaluation/Binary_2_MIX_prune_PRS_AUC_diff_snplists.csv"))

AUC_table_Binary_3 = AUC_table_Binary_3[,c("GWAS_type","approx","trait","target_pop","snplist","weight_AUC")]
AUC_table_Binary_2 = AUC_table_Binary_2[,c("GWAS_type","approx","trait","target_pop","snplist","weight_AUC")]

AUC_table = rbind(AUC_table_Binary_3,AUC_table_Binary_2)

AUC_table$snplist = paste0("snplist_",AUC_table$snplist)
AUC_table$approx = ifelse(AUC_table$approx == "TRUE","Identity","Reference LD")
AUC_table$target_pop = factor(AUC_table$target_pop, levels = c("EUR","EAS","AFR","SAS","AMR"))

AUC_table = AUC_table[,c("target_pop","trait","snplist","weight_AUC")]
colnames(AUC_table) = c("Pop","Trait","Strategy","Metric")

## organize the final table
long_table <- rbind(r2_table,AUC_table)
long_table$Strategy <- factor(long_table$Strategy, levels = c("snplist_1","snplist_2","snplist_3","snplist_4"))
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

write_xlsx(reshaped_table, "TableS11.xlsx")