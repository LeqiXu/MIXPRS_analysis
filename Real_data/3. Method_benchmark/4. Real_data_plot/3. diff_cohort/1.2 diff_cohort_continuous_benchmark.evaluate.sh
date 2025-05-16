## auto and tuning method comparison
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

setwd("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/plot")

# Fixed plot value
all_method_color = c("MIXPRS" = "#6A0DAD","JointPRS" = "#B0003C",
"JointPRS-auto" = "#E65475", "XPASS" = "#FF8C00", "SDPRX" = "#fec44f","PRS-CSx-auto" = "#FDCAC7", 
"PRS-CSx" = "#006400","MUSSEL" = "#8FBC8F","PROSPER" = "#5E92F3","BridgePRS" = "#89CFF0")

my_theme <- theme(
    plot.title = element_text(size=16, face = "bold"),
    text = element_text(size=16),
    axis.text.x = element_text(size=14),
    axis.text.y = element_text(size=14),
    axis.title.x = element_text(size=16, face = "bold"),
    axis.title.y = element_text(size=16),
    strip.text = element_text(face = "bold", size = 14),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey", size = 0.5),
    panel.grid.minor.y = element_line(color = "lightgrey", size = 0.25)
  )

# Continuous plot
# Continuous plot
## Data preprocessing
prs_table = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/evaluation/PRS_update_diff_cohort_r2.csv"))
prs_table = prs_table[,c("pop", "trait", "MIXPRS", "JointPRS_tune_max","SDPRX_auto_2","XPASS_auto_2","PRScsx_tune_max","PROSPER_tune_max","MUSSEL_tune_max","BridgePRS_tune_2")]
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

## Calcualte the relative change of other methods over JointPRS
MIXPRS_ref <- long_table[long_table$method == "MIXPRS", c("pop", "trait", "mean_r2")]
colnames(MIXPRS_ref)[3] <- "MIXPRS_mean_r2"
long_table_with_MIXPRS  <- merge(long_table, MIXPRS_ref, by = c("pop", "trait"))

long_table_with_MIXPRS$relative_change <- with(long_table_with_MIXPRS , 
                                                 (MIXPRS_mean_r2 - mean_r2) / mean_r2)
long_table_with_MIXPRS$relative_change <- long_table_with_MIXPRS$relative_change*100
long_table_with_MIXPRS <- long_table_with_MIXPRS [long_table_with_MIXPRS $method != "MIXPRS",]
long_table_with_MIXPRS <- long_table_with_MIXPRS [,c("pop", "trait","method","relative_change")]
long_table_with_MIXPRS$method = factor(long_table_with_MIXPRS$method, levels = c("JointPRS","XPASS","SDPRX","PRS-CSx","MUSSEL","PROSPER","BridgePRS"))

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