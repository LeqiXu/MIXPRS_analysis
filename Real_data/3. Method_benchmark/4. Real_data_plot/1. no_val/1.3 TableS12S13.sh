## MIXPRS comparison
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(writexl)

setwd("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/table")

## Data preprocessing
GLGC_prs_table = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/evaluation/GLGC_PRS_r2.csv"))
PAGE_prs_table = fread("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/evaluation/PAGE_PRS_r2.csv")
BBJ_prs_table = fread("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/evaluation/BBJ_PRS_r2.csv")

continuous_table = rbind(GLGC_prs_table,PAGE_prs_table)
continuous_table = rbind(continuous_table,BBJ_prs_table)

continuous_table = continuous_table[,c("pop", "trait", "MIX", "JointPRS_auto_max","PRScsx_auto_max","SDPRX_auto_2","XPASS_auto_2")]
colnames(continuous_table) = c("Pop", "Trait", "MIXPRS", "JointPRS-auto","PRS-CSx-auto","SDPRX","XPASS")
continuous_table$trait_type = "Continuous"

# Binary trait
Binary_3_prs_table = fread("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/evaluation/Binary_3_PRS_AUC.csv")
Binary_2_prs_table = fread("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/evaluation/Binary_2_PRS_AUC.csv")

binary_table = rbind(Binary_3_prs_table,Binary_2_prs_table)

binary_table = binary_table[,c("pop", "trait", "MIX", "JointPRS_auto_max","PRScsx_auto_max","SDPRX_auto_2","XPASS_auto_2")]
colnames(binary_table) = c("Pop", "Trait", "MIXPRS", "JointPRS-auto","PRS-CSx-auto","SDPRX","XPASS")
binary_table$trait_type = "Binary"

## All traits
all_table <- rbind(continuous_table,binary_table)
all_table$Pop <- factor(all_table$Pop, levels = c("EUR","EAS","AFR","SAS","AMR"))
all_table$Trait <- factor(all_table$Trait, levels = c("HDL","LDL","TC","logTG","Height","BMI","SBP","DBP","PLT","WBC","NEU","LYM","MON","EOS","RBC","HB","HCT","MCH","MCV","ALT","ALP","GGT","T2D","BrC","CAD","LuC"))

# Round Metric to three decimals and reshape
reshaped_table <- all_table %>%
  select(Pop, Trait, MIXPRS,`JointPRS-auto`,`PRS-CSx-auto`,SDPRX,XPASS) %>%
  arrange(Pop,Trait)

# Check results
print(reshaped_table)

write_xlsx(reshaped_table, "TableS12S13.xlsx")