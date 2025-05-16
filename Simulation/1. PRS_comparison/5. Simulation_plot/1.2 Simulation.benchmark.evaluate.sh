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

setwd("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/plot")

# Continuous plot
# Continuous plot
## Data preprocessing
prs_table = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/evaluation/sim_PRS_real_tune_benchmark_r2.csv"))
prs_table = prs_table[,c("n","pop","p","rhog","sample2","MIXPRS_auto_5","SDPRX_auto_2","XPASS_auto_2","JointPRS_tune_5","PRScsx_tune_5","PROSPER_tune_5","MUSSEL_tune_5","BridgePRS_tune_2")]
colnames(prs_table) = c("n","pop","p","rhog","sample2","MIXPRS","SDPRX","XPASS","JointPRS","PRS-CSx","PROSPER","MUSSEL","BridgePRS")

## Reshape the data to a long format and estimate the mean and sd
long_table <- melt(prs_table, id.vars = c("n","pop","p","rhog","sample2"),
                   variable.name = "method", value.name = "r2")

long_table$p <- factor(long_table$p,
                       levels = c(0.1, 0.01, 0.001, 5e-04),
                       labels = c("p = 0.1", "p = 0.01", "p = 0.001", "p = 5 × 10⁻⁴"))
long_table$pop <- factor(long_table$pop, levels = c("EAS", "AFR", "SAS", "AMR"))
long_table$method = factor(long_table$method, levels = c("MIXPRS","JointPRS","XPASS","SDPRX","PRS-CSx","MUSSEL","PROSPER","BridgePRS"))

## Calcualte the relative change of other methods over JointPRS
long_table <- long_table %>%
  group_by(pop,p,rhog,sample2,method) %>%
  summarise(
    mean_r2 = mean(r2, na.rm = TRUE),
    .groups = 'drop' # This argument drops the grouping structure afterwards
  )

# Calculate improvements relative to MIXPRS
improvement_summary <- long_table %>%
  pivot_wider(names_from = method, values_from = mean_r2) %>%
  pivot_longer(cols = c("JointPRS","XPASS","SDPRX","PRS-CSx","MUSSEL","PROSPER","BridgePRS"),
               names_to = "Compared_Method",
               values_to = "Compared_Metric") %>%
  mutate(Improvement = (MIXPRS - Compared_Metric) / Compared_Metric) %>%
  group_by(pop,sample2,Compared_Method) %>%
  summarize(
    avg_improvement = mean(Improvement, na.rm = TRUE) * 100,
    count = sum(!is.na(Improvement)),
    .groups = "drop"
  )

improvement_summary$avg_improvement <- round(improvement_summary$avg_improvement,2)
improvement_summary$Compared_Method <- factor(improvement_summary$Compared_Method, levels = c("JointPRS","XPASS","SDPRX","PRS-CSx","MUSSEL","PROSPER","BridgePRS"))
improvement_summary <- improvement_summary %>% arrange(pop,sample2,Compared_Method)
# Display the summarized improvement table
print(improvement_summary)
