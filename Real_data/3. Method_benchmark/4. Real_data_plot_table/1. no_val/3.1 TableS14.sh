## MIXPRS comparison
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(writexl)

setwd("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/table")

# Continuous trait
GLGC_prs_table = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/evaluation/GLGC_PUMAS_EN_PRS_r2.csv"))
colnames(GLGC_prs_table) = c("Pop","Trait","MIXPRS","PUMAS-EN","PUMAS-EN_paper")
GLGC_prs_table$Pop <- factor(GLGC_prs_table$Pop, levels = c("EUR","EAS","AFR","SAS","AMR"))
GLGC_prs_table$Trait <- factor(GLGC_prs_table$Trait, levels = c("HDL","LDL","TC","logTG","Height","BMI","SBP","DBP","PLT","WBC","NEU","LYM","MON","EOS","RBC","HB","HCT","MCH","MCV","ALT","ALP","GGT","T2D","BrC","CAD","LuC"))

reshaped_table <- GLGC_prs_table %>%
  arrange(Pop,Trait)

# Check results
print(reshaped_table)

write_xlsx(reshaped_table, "TableS14.xlsx")