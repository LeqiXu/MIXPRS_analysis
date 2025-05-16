library(data.table)
library(dplyr)
library(tidyr)
library(writexl)

setwd("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/table")

## Data preprocessing
prs_table = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/evaluation/sim_PRS_real_MIXPRS_benefit_r2.csv"))
long_table <- melt(prs_table, id.vars = c("n", "pop", "p", "rhog", "sample1", "sample2"),
                   variable.name = "method", value.name = "r2")
long_table$pop <- factor(long_table$pop, levels = c("EUR","EAS","AFR","SAS","AMR"))
long_table$p <- factor(long_table$p,
                       levels = c(0.1, 0.01, 0.001, 5e-04),
                       labels = c("p = 0.1", "p = 0.01", "p = 0.001", "p = 5 × 10⁻⁴"))
long_table$method <- dplyr::recode(
  long_table$method,
  "JointPRS_auto"  = "JointPRS-auto",
  "SDPRX" = "SDPRX",
  "MIXPRS_prune_TRUE"  = "MIXPRS"
)

long_table$pop <- factor(long_table$pop, levels = c("EAS", "AFR", "SAS", "AMR"))
long_table$method = factor(long_table$method, levels = c("MIXPRS","JointPRS-auto","SDPRX"))
long_table = long_table[,c("n","pop","p","sample2","method","r2")]

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

write_xlsx(reshaped_table, "TableS6.xlsx")