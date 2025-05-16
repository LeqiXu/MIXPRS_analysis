## MIXPRS comparison
library(data.table)
library(dplyr)
library(tidyr)
library(writexl)

setwd("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/table")

corr_table = fread("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/residual/simulation_residual_corr.txt")
corr_table$pop = factor(corr_table$pop, levels = c("EUR","EAS","AFR","SAS","AMR"))
corr_table$type <- dplyr::recode(
  corr_table$type,
  "full_snplist"  = "RefLD_Full",
  "prune_snplist" = "RefLD_Prune",
  "prune_snplist_ind_approx"  = "Identity_Prune"
)
corr_table$type = factor(corr_table$type, levels = c("RefLD_Full","RefLD_Prune","Identity_Prune"))
corr_table$p <- factor(corr_table$p,
                       levels = c(0.1, 0.01, 0.001, 5e-04),
                       labels = c("p = 0.1", "p = 0.01", "p = 0.001", "p = 5 × 10⁻⁴"))

summary_corr_table<- corr_table %>%
  group_by(pop, p, type) %>%
  summarise(mean_corr = mean(corr), .groups = "drop")

# Structure the data to have separate columns for mean_r2 for each combination of population and sample2:
reshaped_table <- summary_corr_table %>%
  select(p, type, pop, mean_corr) %>%
  pivot_wider(names_from = pop, values_from = mean_corr) %>%
  arrange(p, type)

# Optional: Replace NA with zeros if needed
reshaped_table[is.na(reshaped_table)] <- 0

# Preview the reshaped table
print(reshaped_table)

write_xlsx(reshaped_table, "TableS7.xlsx")