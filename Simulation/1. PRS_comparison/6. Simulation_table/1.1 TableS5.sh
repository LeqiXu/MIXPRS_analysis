library(data.table)
library(dplyr)
library(tidyr)
library(writexl)

setwd("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/table")

## Data preprocessing
prs_table = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/evaluation/sim_PRS_real_tune_benchmark_r2.csv"))
prs_table = prs_table[,c("n","pop","p","sample2","MIXPRS_auto_5","SDPRX_auto_2","XPASS_auto_2","JointPRS_tune_5","PRScsx_tune_5","PROSPER_tune_5","MUSSEL_tune_5","BridgePRS_tune_2")]
colnames(prs_table) = c("n","pop","p","sample2","MIXPRS","SDPRX","XPASS","JointPRS","PRS-CSx","PROSPER","MUSSEL","BridgePRS")

## Reshape the data to a long format and estimate the mean and sd
long_table <- melt(prs_table, id.vars = c("n","pop","p","sample2"),
                   variable.name = "method", value.name = "r2")

long_table$p <- factor(long_table$p,
                       levels = c(0.1, 0.01, 0.001, 5e-04),
                       labels = c("p = 0.1", "p = 0.01", "p = 0.001", "p = 5 × 10⁻⁴"))
long_table$pop <- factor(long_table$pop, levels = c("EAS", "AFR", "SAS", "AMR"))
long_table$method = factor(long_table$method, levels = c("MIXPRS","JointPRS","XPASS","SDPRX","PRS-CSx","MUSSEL","PROSPER","BridgePRS"))

## Calcualte the relative change of other methods over JointPRS
long_table <- long_table %>%
  group_by(pop,p,sample2,method) %>%
  summarise(
    mean_r2 = mean(r2, na.rm = TRUE),
    .groups = 'drop' # This argument drops the grouping structure afterwards
  )

# Structure the data to have separate columns for mean_r2 for each combination of population and sample2:
reshaped_table <- long_table %>%
  unite("pop_sample", pop, sample2, sep = "_") %>% # combines population and sample2 into one column
  select(p, method, pop_sample, mean_r2) %>%
  pivot_wider(names_from = pop_sample, values_from = mean_r2) %>%
  arrange(p, method)

# Optional: Replace NA with zeros if needed
reshaped_table[is.na(reshaped_table)] <- 0

# Preview the reshaped table
print(reshaped_table)

write_xlsx(reshaped_table, "TableS5.xlsx")