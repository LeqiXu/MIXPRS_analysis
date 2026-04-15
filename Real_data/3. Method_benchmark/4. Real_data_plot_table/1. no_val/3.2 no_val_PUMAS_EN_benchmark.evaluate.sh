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
library(ggpattern)

# Continuous trait
GLGC_prs_table = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/evaluation/GLGC_PUMAS_EN_PRS_r2.csv"))
colnames(GLGC_prs_table) = c("pop","trait","MIXPRS","PUMAS-EN","PUMAS-EN_paper")

## Reshape the data to a long format and estimate the mean and sd
long_table <- melt(GLGC_prs_table, id.vars = c("pop", "trait"),
                   variable.name = "method", value.name = "metric")
long_table$trait = factor(long_table$trait, levels = c("HDL","LDL","TC","logTG"))
long_table$pop <- factor(long_table$pop, levels = c("EAS","AFR","SAS","AMR"))
long_table$method = factor(long_table$method, levels = c("MIXPRS","PUMAS-EN","PUMAS-EN_paper"))

# Calculate improvements relative to MIXPRS
improvement_summary <- long_table %>%
  pivot_wider(names_from = method, values_from = metric) %>%
  pivot_longer(cols = c("PUMAS-EN","PUMAS-EN_paper"),
               names_to = "Compared_Method",
               values_to = "Compared_Metric") %>%
  mutate(Improvement = (MIXPRS - Compared_Metric) / Compared_Metric) %>%
  group_by(pop, Compared_Method) %>%
  summarize(
    avg_improvement = mean(Improvement, na.rm = TRUE) * 100,
    count = sum(!is.na(Improvement)),
    .groups = "drop"
  )

improvement_summary$avg_improvement <- round(improvement_summary$avg_improvement,2)
# Display the summarized improvement table
print(improvement_summary)
