# Step 2: Organize beta by chr pop for each param in each scenario
library(data.table)

# --- Define Paths ---
# This is the main output directory you specified in your SLURM script
base_result_dir <- "/home/yd357/pi_paths/pi_zhao/MIXPRS_Revise/Project1_15K_80K/result/add_benchmark/25.10.4/sim_data/SDPRX"
# We will create a new subdirectory to store the combined beta files
output_dir <- file.path(base_result_dir, "beta_organized")

if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}

# --- Define Parameters ---
rhog <- "0.8"
pop1 <- "EUR"
sample_sizes_2 <- c("15K", "80K") # Your loop for sample_size_2
sim_is <- 1:5
ps <- c("0.1", "0.01", "0.001", "5e-04")
pop2s <- c("EAS", "AFR", "SAS", "AMR") # Your loop for pop2

# --- Main Loop ---
# Loop over all simulation parameters
for (sample_size_2 in sample_sizes_2) {
    for (sim_i in sim_is) {
        for (p in ps) {
            for (pop2 in pop2s) {

                cat(paste("Processing: Sim", sim_i, "- p", p, "- Sample", sample_size_2, "- Pop", pop2, "\n"))

                # --- 1. Process Betas for Pop2 (from _2.txt files) ---
                SDPRX_all_pop2 <- data.table()
                all_pop2_files_found <- TRUE
                
                for (chr in 1:22) {
                    # Construct the input file name based on your job_run_sdprx.sh OUT_NAME
                    # sim${SIM_I}_p${P}_rho${RHOG}_${SAMPLE_SIZE_2}_${POP1}_${POP2}_SDPRX_chr${CHR}_2.txt
                    input_file_2 <- file.path(base_result_dir,
                                              paste0("sim", sim_i, "_p", p, "_rho", rhog, "_", sample_size_2,
                                                     "_", pop1, "_", pop2, "_SDPRX_chr", chr, "_2.txt"))
                    
                    if (!file.exists(input_file_2)) {
                        cat(paste("  WARNING: Missing file:", input_file_2, "\n"))
                        all_pop2_files_found <- FALSE
                        next
                    }
                    
                    SDPRX_pop_chr <- fread(input_file_2)
                    # The file contains betas for pop2
                    names(SDPRX_pop_chr) <- c("rsID", "A1", pop2)
                    
                    SDPRX_all_pop2 <- rbind(SDPRX_all_pop2, SDPRX_pop_chr)
                }
                
                # Write the combined file for pop2
                if (nrow(SDPRX_all_pop2) > 0) {
                    output_file_pop2 <- file.path(output_dir,
                                                  paste0("sim", sim_i, "_p", p, "_rho", rhog, "_", sample_size_2,
                                                         "_", pop1, "_", pop2, "_beta_", pop2, ".txt"))
                    write.table(SDPRX_all_pop2, output_file_pop2, quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)
                } else if (all_pop2_files_found) {
                     cat(paste("  WARNING: No data found for pop2", pop2, "combination.\n"))
                }

                # --- 2. Process Betas for Pop1 (EUR) (from _1.txt files) ---
                SDPRX_all_eur <- data.table()
                all_eur_files_found <- TRUE
                
                for (chr in 1:22) {
                    # Construct the input file name based on your job_run_sdprx.sh OUT_NAME
                    # sim${SIM_I}_p${P}_rho${RHOG}_${SAMPLE_SIZE_2}_${POP1}_${POP2}_SDPRX_chr${CHR}_1.txt
                    input_file_1 <- file.path(base_result_dir,
                                              paste0("sim", sim_i, "_p", p, "_rho", rhog, "_", sample_size_2,
                                                     "_", pop1, "_", pop2, "_SDPRX_chr", chr, "_1.txt"))
                    
                    if (!file.exists(input_file_1)) {
                        cat(paste("  WARNING: Missing file:", input_file_1, "\n"))
                        all_eur_files_found <- FALSE
                        next
                    }
                    
                    SDPRX_pop_chr <- fread(input_file_1)
                    # The file contains betas for pop1 (EUR)
                    names(SDPRX_pop_chr) <- c("rsID", "A1", pop1)
                    
                    SDPRX_all_eur <- rbind(SDPRX_all_eur, SDPRX_pop_chr)
                }
                
                # Write the combined file for EUR
                if (nrow(SDPRX_all_eur) > 0) {
                    output_file_eur <- file.path(output_dir,
                                                 paste0("sim", sim_i, "_p", p, "_rho", rhog, "_", sample_size_2,
                                                        "_", pop1, "_", pop2, "_beta_", pop1, ".txt"))
                    write.table(SDPRX_all_eur, output_file_eur, quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)
                } else if (all_eur_files_found) {
                     cat(paste("  WARNING: No data found for pop1 (EUR) in", pop2, "combination.\n"))
                }
            }
        }
    }
}

cat("--- Result organization complete. ---\n")
