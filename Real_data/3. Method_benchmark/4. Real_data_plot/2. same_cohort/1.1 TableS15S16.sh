## MIXPRS comparison
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(writexl)

setwd("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/table")

## Data preprocessing
prs_table1 = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/evaluation/PRS_update_same_cohort_r2.csv"))
prs_table2 = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/evaluation/PRS_update_same_cohort_auc.csv"))
prs_table = rbind(prs_table1,prs_table2)
prs_table = prs_table[,c("pop", "trait", "MIX" , "JointPRS_tune_max","SDPRX_auto_2","XPASS_auto_2","PRScsx_tune_max","PROSPER_tune_max","MUSSEL_tune_max","BridgePRS_tune_2")]
colnames(prs_table) = c("Pop", "Trait", "MIXPRS", "JointPRS","SDPRX","XPASS","PRS-CSx","PROSPER","MUSSEL","BridgePRS")

## Reshape the data to a long format and estimate the mean and sd
prs_mean_table <- prs_table %>%
  group_by(Pop, Trait) %>%
  summarise(
    MIXPRS = mean(MIXPRS, na.rm = TRUE),
    JointPRS = mean(JointPRS, na.rm = TRUE),
    SDPRX = mean(SDPRX, na.rm = TRUE),
    XPASS = mean(XPASS, na.rm = TRUE),
    `PRS-CSx` = mean(`PRS-CSx`, na.rm = TRUE),
    PROSPER = mean(PROSPER, na.rm = TRUE),
    MUSSEL = mean(MUSSEL, na.rm = TRUE),
    BridgePRS = mean(BridgePRS, na.rm = TRUE),
    .groups = 'drop' # This argument drops the grouping structure afterwards
  )

prs_mean_table$Pop <- factor(prs_mean_table$Pop, levels = c("EUR","EAS","AFR","SAS","AMR"))
prs_mean_table$Trait <- factor(prs_mean_table$Trait, levels = c("HDL","LDL","TC","logTG","Height","BMI","SBP","DBP","PLT","WBC","NEU","LYM","MON","EOS","RBC","HB","HCT","MCH","MCV","ALT","ALP","GGT","T2D","BrC","CAD","LuC"))

reshaped_table <- prs_mean_table %>%
  arrange(Pop,Trait)

# Check results
print(reshaped_table)

write_xlsx(reshaped_table, "TableS15S16.xlsx")