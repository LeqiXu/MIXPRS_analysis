## MIXPRS comparison
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(writexl)

setwd("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/table")

## Data preprocessing
prs_table = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/evaluation/PRS_update_diff_cohort_r2.csv"))
prs_table = prs_table[,c("pop", "trait", "MIXPRS", "JointPRS_tune_max","SDPRX_auto_2","XPASS_auto_2","PRScsx_tune_max","PROSPER_tune_max","MUSSEL_tune_max","BridgePRS_tune_2")]
colnames(prs_table) = c("Pop", "Trait", "MIXPRS", "JointPRS","SDPRX","XPASS","PRS-CSx","PROSPER","MUSSEL","BridgePRS")

prs_table$Pop <- factor(prs_table$Pop, levels = c("EUR","EAS","AFR","SAS","AMR"))
prs_table$Trait <- factor(prs_table$Trait, levels = c("HDL","LDL","TC","logTG","Height","BMI","SBP","DBP","PLT","WBC","NEU","LYM","MON","EOS","RBC","HB","HCT","MCH","MCV","ALT","ALP","GGT","T2D","BrC","CAD","LuC"))

reshaped_table <- prs_table %>%
  arrange(Pop,Trait)

# Check results
print(reshaped_table)

write_xlsx(reshaped_table, "TableS17.xlsx")