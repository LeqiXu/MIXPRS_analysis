# ---- Setup ----
library(data.table)
library(stringr)
library(ggplot2)
library(ggpubr)

setwd("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/plot")

my_theme <- theme(
  # shrink axis titles
  axis.title.x = element_text(size = 6, face = "bold"),
  axis.title.y = element_text(size = 6, face = "bold"),
  # shrink tick labels
  axis.text.x  = element_text(size = 4, face = "bold"),
  axis.text.y  = element_text(size = 4, face = "bold"),
  # leave all other text at your preferred size
  plot.title   = element_text(size = 8, face = "bold"),
  strip.text   = element_text(size = 6, face = "bold"),
  panel.grid.major.x = element_blank(),
  panel.grid.minor.x = element_blank(),
  panel.grid.major.y = element_blank(),
  panel.grid.minor.y = element_blank(),
  legend.position  = "top",
  legend.title     = element_text(size = 7, face = "bold"),     # title size
  legend.text      = element_text(size = 5)
)

# Root directory and parameters
base_dir <- "/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/residual"

h2   <- "0.4"
rhog <- "0.8" 
nlab <- "15K"

pops     <- c("EAS", "AFR", "SAS", "AMR")
sims     <- c(1:5)
p_vals   <- c(0.1, 0.01, 0.001, 5e-04)


# Filename pattern:
# Pattern below: {POP}_sim{SIM}_h2{h2}_p{P}_rhog{rhog}_{nlab}_{POP}_sim_GWAS_residual.txt
file_for <- function(pop, sim, p) {
  p_str <- if (p == 5e-04) "5e-04" else as.character(p)
  fname <- sprintf("%s_sim%d_h2%s_p%s_rhog%s_%s_%s_sim_GWAS_residual.txt",
                   pop, sim, h2, p_str, rhog, nlab, pop)
  file.path(base_dir, fname)
}

# ---- Helper: compute quantile QQ summary for one numeric vector ----
qq_summary <- function(x, n = 999, trim = 0.01) {
  # Standardize to focus on Gaussian shape (location/scale invariant)
  x <- as.numeric(x)
  x <- (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)

  # Use interior probs to avoid extreme tails instability on huge N
  probs <- seq(trim, 1 - trim, length.out = n)
  data.table(
    prob = probs,
    sample_q = quantile(x, probs, na.rm = TRUE, type = 8),
    theo_q   = qnorm(probs)  # standard normal
  )
}

# ---- Ingest and summarize all files ----
all_summaries <- rbindlist(
  lapply(pops, function(pop) {
    rbindlist(lapply(p_vals, function(p) {
      rbindlist(lapply(sims, function(sim) {
        f <- file_for(pop, sim, p)
        if (!file.exists(f)) {
          warning("Missing file: ", f)
          return(NULL)
        }
        dt <- tryCatch(fread(f, select = "Residual"), error = function(e) NULL)
        if (is.null(dt) || !"Residual" %in% names(dt)) return(NULL)

        qq <- qq_summary(dt$Residual)
        qq[, `:=`(pop = pop, p = p, sim = sim)]
        qq
      }), use.names = TRUE, fill = TRUE)
    }), use.names = TRUE, fill = TRUE)
  }), use.names = TRUE, fill = TRUE
)

# If nothing was read (path mismatch), stop with a helpful message
if (nrow(all_summaries) == 0) {
  stop("No summaries computed. Check base_dir and filename pattern in file_for().")
}

# ---- Aggregate across simulations: median line + 10–90% ribbon ----
agg <- all_summaries[
  , .(
      q_median = median(sample_q, na.rm = TRUE),
      q_lo10   = quantile(sample_q, 0.10, na.rm = TRUE, type = 8),
      q_hi90   = quantile(sample_q, 0.90, na.rm = TRUE, type = 8)
    ),
    by = .(pop, p, prob, theo_q)
]

# Make p a factor for nice facet labels
agg[, p_fac := factor(p, levels = p_vals, labels = paste0("p=", c("0.1","0.01","0.001","5e-04")))]
agg[, pop := factor(pop, levels = pops)]

# ---- Plot: one clean figure ----
p <- ggplot(agg, aes(x = theo_q)) +
  geom_ribbon(aes(ymin = q_lo10, ymax = q_hi90), alpha = 0.18) +
  geom_line(aes(y = q_median), linewidth = 0.8) +
  geom_abline(slope = 1, intercept = 0, color = "red", linewidth = 0.3) +
  facet_grid(p_fac ~ pop, labeller = label_value) +
  labs(
    x = "Theoretical quantile for N(0,1)",
    y = "Sample quantile for standardized residuals"
  ) +
  theme_classic() +
  my_theme

print(p)

# Save figure (change path/size as needed)
ggsave(filename = "FigureS2.pdf",
       plot = p,
       width = 6.5, 
       height = 6,
       units = "in",
       dpi = 300,
       device   = cairo_pdf) 