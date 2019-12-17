# ========================================================
# Runner script for simulations and analysis
# ========================================================

##########################################################
# Set up working environment
project_dir <- "~/Github/RTG-2019-2H"
setwd(project_dir)
source("R/setup.R", local=FALSE, echo=TRUE)
source("R/simulation.R", local=FALSE, echo=TRUE)
source("R/estimator.R", local=FALSE, echo=TRUE)
source("R/reporting.R", local=FALSE, echo=TRUE)
source("R/result_export.R", local=FALSE, echo=TRUE)
# source("R/propensity_score.R", local=FALSE, echo=TRUE)
project_dir <- "~/Github/RTG-2019-2H"
setwd(project_dir)
output_file = file(paste0("outputs/R_log_", format(Sys.time(), 
                   "%Y%m%d%H%M"), ".txt"), "w")

##########################################################
# Simulation initialization and loop code

# Load cached datasets and run results
n_methods = 16
n_sim = 500
n_dgps = 6
# n_columns = 15
n_df = n_sim*n_dgps*n_methods

# Create outputs file to store results
# output = matrix(rep(-Inf, n_sim*n_dgps*n_columns*n_methods),
#                 ncol = n_columns)
output_df <- data.frame(sim_num=rep(NA, n_df), 
                        estimator=rep("", n_df), 
                        method=rep("", n_df), 
                        outcome=rep("", n_df), 
                        effect_type=rep("", n_df), 
                        unlabeled_mechanism=rep("", n_df), 
                        pct_unlab=rep("", n_df), 
                        confoundedness=rep("", n_df), 
                        other_notes=rep("", n_df), 
                        sample_size=rep(NA, n_df), 
                        true_tau=rep(NA, n_df), 
                        ATE=rep(NA, n_df), 
                        ATE_CI_lb=rep(NA, n_df), 
                        ATE_CI_ub=rep(NA, n_df), 
                        rmse=rep(NA, n_df), 
                        cover=rep(NA, n_df), 
                        stringsAsFactors=FALSE)

for (i in 1:n_sim){
  t = proc.time()[3]
  
  ##########################################################
  # First DGP
  dgp_num = 1; startrow = 1 + (i-1)*n_dgps*n_methods; endrow = startrow + n_methods - 1
  sample = 5000; lin = "lin"; hom = "hom"
  missing_type = "MCAR"; missing_pct = "90"; conf = "unconf"
  subfolder_name = paste(lin, hom, sample, missing_type, missing_pct, conf, sep="_")
  file_name = paste0("df", i, ".csv")
  data_file <- file.path(project_dir, "data", "simulations", subfolder_name, file_name)
  df <- read.csv(data_file)
  
  # Run all estimates on the DGP loaded into memory
  output_df[startrow:endrow, ] <- simulation_estimates(
    i, df, sample, lin, hom, missing_type, missing_pct, conf
  )
  
  ##########################################################
  # Second DGP
  dgp_num = 2; startrow = 1 + endrow; endrow = endrow + n_methods
  sample = 5000; lin = "lin"; hom = "hom"
  missing_type = "MCAR"; missing_pct = "90"; conf = "conf"
  subfolder_name = paste(lin, hom, sample, missing_type, missing_pct, conf, sep="_")
  file_name = paste0("df", i, ".csv")
  data_file <- file.path(project_dir, "data", "simulations", subfolder_name, file_name)
  df <- read.csv(data_file)
  
  # Run all estimates on the DGP loaded into memory
  output_df[startrow:endrow, ] <- simulation_estimates(
    i, df, sample, lin, hom, missing_type, missing_pct, conf
  )
  
  ##########################################################
  # Third DGP
  dgp_num = 3; startrow = 1 + endrow; endrow = endrow + n_methods
  sample = 1000; lin = "lin"; hom = "hom"
  missing_type = "MCAR"; missing_pct = "90"; conf = "unconf"
  subfolder_name = paste(lin, hom, sample, missing_type, missing_pct, conf, sep="_")
  file_name = paste0("df", i, ".csv")
  data_file <- file.path(project_dir, "data", "simulations", subfolder_name, file_name)
  df <- read.csv(data_file)
  
  # Run all estimates on the DGP loaded into memory
  output_df[startrow:endrow, ] <- simulation_estimates(
    i, df, sample, lin, hom, missing_type, missing_pct, conf
  )
  
  ##########################################################
  # Fourth DGP
  dgp_num = 4; startrow = 1 + endrow; endrow = endrow + n_methods
  sample = 1000; lin = "lin"; hom = "hom"
  missing_type = "MCAR"; missing_pct = "90"; conf = "conf"
  subfolder_name = paste(lin, hom, sample, missing_type, missing_pct, conf, sep="_")
  file_name = paste0("df", i, ".csv")
  data_file <- file.path(project_dir, "data", "simulations", subfolder_name, file_name)
  df <- read.csv(data_file)
  
  # Run all estimates on the DGP loaded into memory
  output_df[startrow:endrow, ] <- simulation_estimates(
    i, df, sample, lin, hom, missing_type, missing_pct, conf
  )
  
  ##########################################################
  # Fifth DGP
  dgp_num = 5; startrow = 1 + endrow; endrow = endrow + n_methods
  sample = 500; lin = "lin"; hom = "hom"
  missing_type = "MCAR"; missing_pct = "90"; conf = "unconf"
  subfolder_name = paste(lin, hom, sample, missing_type, missing_pct, conf, sep="_")
  file_name = paste0("df", i, ".csv")
  data_file <- file.path(project_dir, "data", "simulations", subfolder_name, file_name)
  df <- read.csv(data_file)
  
  # Run all estimates on the DGP loaded into memory
  output_df[startrow:endrow, ] <- simulation_estimates(
    i, df, sample, lin, hom, missing_type, missing_pct, conf
  )
  
  ##########################################################
  # Sixth DGP
  dgp_num = 6; startrow = 1 + endrow; endrow = endrow + n_methods
  sample = 500; lin = "lin"; hom = "hom"
  missing_type = "MCAR"; missing_pct = "90"; conf = "conf"
  subfolder_name = paste(lin, hom, sample, missing_type, missing_pct, conf, sep="_")
  file_name = paste0("df", i, ".csv")
  data_file <- file.path(project_dir, "data", "simulations", subfolder_name, file_name)
  df <- read.csv(data_file)
  
  # Run all estimates on the DGP loaded into memory
  output_df[startrow:endrow, ] <- simulation_estimates(
    i, df, sample, lin, hom, missing_type, missing_pct, conf
  )

  t = proc.time()[3] - t
  cat("==================================================\n", file = output_file)
  cat(paste("Run time:", t, "\n"), file = output_file)
  cat(paste("Sim", i, "of", n_sim, "complete\n\n"), file = output_file)
  # cat("Sim", i, "of", n_sim, "complete\n")
}
close(output_file)

##########################################################
# Save simulation output

# Write overall simulation table for later access
save_simulation_output(output_df, "sim_")

# Collapse results into table of averages
mean_simulations <- summarize_simulation_results(output_df)

# Format table outputs
mean_simulations <- format_simulation_outputs(mean_simulations)

# Write overall simulation summary to a CSV file
save_simulation_output(mean_simulations, "summary_table_")

##########################################################
# Save as Latex tables for use in paper / presentation

### First output table
#  - Semi-supervised
#  - MCAR 90%
#  - IPW with logistic regression
#  - Linear Outcome
#  - Homogeneous Treatment Effect
xtable_df <- prepare_per_estimator_table(mean_simulations, "IPW")
xtable_output <- format_per_estimator_xtable(xtable_df, n_sim, "MCAR")
save_per_estimator_to_latex(xtable_output, "IPW")

### Second output table
#  - Semi-supervised
#  - MCAR 90%
#  - IPW with BART
#  - Linear Outcome
#  - Homogeneous Treatment Effect
xtable_df <- prepare_per_estimator_table(mean_simulations, "IPW_BART")
xtable_output <- format_per_estimator_xtable(xtable_df, n_sim, "MCAR")
save_per_estimator_to_latex(xtable_output, "IPW_BART")

### Third output table
#  - Semi-supervised
#  - MCAR 90%
#  - IPW with BART
#  - Linear Outcome
#  - Homogeneous Treatment Effect
xtable_df <- prepare_per_estimator_table(mean_simulations, "TMLE")
xtable_output <- format_per_estimator_xtable(xtable_df, n_sim, "MCAR")
save_per_estimator_to_latex(xtable_output, "TMLE")

### Fourth output table
#  - Semi-supervised
#  - MCAR 90%
#  - BCF
#  - Linear Outcome
#  - Homogeneous Treatment Effect
xtable_df <- prepare_per_estimator_table(mean_simulations, "BCF")
xtable_output <- format_per_estimator_xtable(xtable_df, n_sim, "MCAR")
save_per_estimator_to_latex(xtable_output, "BCF")

### Fifth output table
#  - Complete case
#  - MCAR 90%
#  - TMLE
#  - Linear Outcome
#  - Homogeneous Treatment Effect
xtable_df <- prepare_all_estimator_table(mean_simulations, "Complete Case")
xtable_output <- format_all_estimator_xtable(xtable_df, n_sim, "MCAR")
save_all_estimator_to_latex(xtable_output, "CC")

### Sixth output table
#  - Semi-supervised
#  - MCAR 90%
#  - TMLE
#  - Linear Outcome
#  - Homogeneous Treatment Effect
xtable_df <- prepare_all_estimator_table(mean_simulations, "Semi-supervised")
xtable_output <- format_all_estimator_xtable(xtable_df, n_sim, "MCAR")
save_all_estimator_to_latex(xtable_output, "SSL")
