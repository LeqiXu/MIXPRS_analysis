## auto and tune method comparison
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

setwd("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/table")

# Continuous plot
## Data preprocessing
prs_table1 = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/evaluation/PRS_update_same_cohort_r2.csv"))
prs_table2 = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/evaluation/PRS_update_same_cohort_auc.csv"))
prs_table = rbind(prs_table1,prs_table2)
colnames(prs_table) = c("pop", "trait", "MIXPRS", "JointPRS","SDPRX","XPASS","PRS-CSx","PROSPER","MUSSEL","BridgePRS")

## Reshape the data to a long format and estimate the mean and sd
long_table <- melt(prs_table, id.vars = c("pop", "trait"),
                   variable.name = "method", value.name = "r2")
long_table$trait = factor(long_table$trait, levels = c("HDL","LDL","TC","logTG","Height","BMI","SBP","DBP","PLT","WBC","NEU","LYM","MON","EOS","RBC","HB","HCT","MCH","MCV","ALT","ALP","GGT"))
long_table$pop <- factor(long_table$pop, levels = c("EAS","AFR","SAS","AMR"))
long_table$method = factor(long_table$method, levels = c("MIXPRS","JointPRS","XPASS","SDPRX","PRS-CSx","MUSSEL","PROSPER","BridgePRS"))
long_table$r2[which(long_table$r2 == 0)] = NA

long_table <- long_table %>%
  group_by(pop, trait, method) %>%
  summarise(
    mean_r2 = mean(r2, na.rm = TRUE),
    sd_r2 = sd(r2, na.rm = TRUE),
    .groups = 'drop' # This argument drops the grouping structure afterwards
  )

GLGC_trait = c("HDL","LDL","TC","logTG")
PAGE_trait = c("Height","BMI","SBP","DBP","PLT")
BBJ_trait = c("WBC","NEU","LYM","MON","EOS","RBC","HB","HCT","MCH","MCV","ALT","ALP","GGT")

## Calcualte the relative change of other methods over MIXPRS
MIXPRS_ref <- long_table[long_table$method == "MIXPRS", c("pop", "trait", "mean_r2")]
colnames(MIXPRS_ref)[3] <- "MIXPRS_mean_r2"
long_table_with_MIXPRS  <- merge(long_table, MIXPRS_ref, by = c("pop", "trait"))

long_table_with_MIXPRS$relative_change <- with(long_table_with_MIXPRS , 
                                                 (MIXPRS_mean_r2 - mean_r2) / mean_r2)
long_table_with_MIXPRS$relative_change <- long_table_with_MIXPRS$relative_change*100
long_table_with_MIXPRS <- long_table_with_MIXPRS [long_table_with_MIXPRS $method != "MIXPRS",]
long_table_with_MIXPRS <- long_table_with_MIXPRS [,c("pop", "trait","method","relative_change")]
long_table_with_MIXPRS$method = factor(long_table_with_MIXPRS$method, levels = c("JointPRS", "XPASS","SDPRX","PRS-CSx","MUSSEL","PROSPER","BridgePRS"))

# Compute relative improvement of MIXPRS over each method across traits by population
relative_improvement_summary <- long_table_with_MIXPRS %>%
  group_by(pop, method) %>%
  summarise(
    mean_relative_improvement = mean(relative_change, na.rm = TRUE),  # negative because original calculation is other - MIXPRS
    num_traits = sum(!is.na(relative_change)),
    .groups = 'drop'
  )
relative_improvement_summary$mean_relative_improvement <- round(relative_improvement_summary$mean_relative_improvement,2)

# Display results clearly:
print(relative_improvement_summary)


# Convert to wide format and add percentage symbol
improvement_table_wide <- relative_improvement_summary %>%
  select(pop, method, mean_relative_improvement) %>%
  pivot_wider(
    names_from = method,
    values_from = mean_relative_improvement
  ) %>%
  mutate(across(-pop, ~ paste0(.x, "%"))) %>%
  arrange(pop)

# Display the formatted wide table
print(improvement_table_wide)

write_xlsx(improvement_table_wide, "Table12.2.xlsx")