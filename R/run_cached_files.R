# ========================================================
# Runner script for cached simulated datasets
# ========================================================

##########################################################
# Set up working environment
require(here)
project_dir <- here()
source(file.path(project_dir, "R", "setup.R"), local=FALSE, echo=TRUE)
source(file.path(project_dir, "R", "simulation.R"), local=FALSE, echo=TRUE)
source(file.path(project_dir, "R", "estimator.R"), local=FALSE, echo=TRUE)
source(file.path(project_dir, "R", "reporting.R"), local=FALSE, echo=TRUE)
source(file.path(project_dir, "R", "result_export.R"), local=FALSE, echo=TRUE)
datestamp = format(Sys.time(), "%Y%m%d%H%M")
logfile_name = file.path(project_dir, "log", paste0("R_batch_log_", datestamp, ".log"))
output_file = file(logfile_name, "w")

##########################################################
# Simulation initialization and loop code

# Load cached datasets and run results
n_methods = 12
n_sim = 500
n_dgps = 6
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
  
  # ##########################################################
  # # Third DGP
  # dgp_num = 3; startrow = 1 + endrow; endrow = endrow + n_methods
  # sample = 1000; lin = "lin"; hom = "hom"
  # missing_type = "MCAR"; missing_pct = "90"; conf = "unconf"
  # subfolder_name = paste(lin, hom, sample, missing_type, missing_pct, conf, sep="_")
  # file_name = paste0("df", i, ".csv")
  # data_file <- file.path(project_dir, "data", "simulations", subfolder_name, file_name)
  # df <- read.csv(data_file)
  # 
  # # Run all estimates on the DGP loaded into memory
  # output_df[startrow:endrow, ] <- simulation_estimates(
  #   i, df, sample, lin, hom, missing_type, missing_pct, conf
  # )
  # 
  # ##########################################################
  # # Fourth DGP
  # dgp_num = 4; startrow = 1 + endrow; endrow = endrow + n_methods
  # sample = 1000; lin = "lin"; hom = "hom"
  # missing_type = "MCAR"; missing_pct = "90"; conf = "conf"
  # subfolder_name = paste(lin, hom, sample, missing_type, missing_pct, conf, sep="_")
  # file_name = paste0("df", i, ".csv")
  # data_file <- file.path(project_dir, "data", "simulations", subfolder_name, file_name)
  # df <- read.csv(data_file)
  # 
  # # Run all estimates on the DGP loaded into memory
  # output_df[startrow:endrow, ] <- simulation_estimates(
  #   i, df, sample, lin, hom, missing_type, missing_pct, conf
  # )
  # 
  # ##########################################################
  # # Fifth DGP
  # dgp_num = 5; startrow = 1 + endrow; endrow = endrow + n_methods
  # sample = 500; lin = "lin"; hom = "hom"
  # missing_type = "MCAR"; missing_pct = "90"; conf = "unconf"
  # subfolder_name = paste(lin, hom, sample, missing_type, missing_pct, conf, sep="_")
  # file_name = paste0("df", i, ".csv")
  # data_file <- file.path(project_dir, "data", "simulations", subfolder_name, file_name)
  # df <- read.csv(data_file)
  # 
  # # Run all estimates on the DGP loaded into memory
  # output_df[startrow:endrow, ] <- simulation_estimates(
  #   i, df, sample, lin, hom, missing_type, missing_pct, conf
  # )
  # 
  # ##########################################################
  # # Sixth DGP
  # dgp_num = 6; startrow = 1 + endrow; endrow = endrow + n_methods
  # sample = 500; lin = "lin"; hom = "hom"
  # missing_type = "MCAR"; missing_pct = "90"; conf = "conf"
  # subfolder_name = paste(lin, hom, sample, missing_type, missing_pct, conf, sep="_")
  # file_name = paste0("df", i, ".csv")
  # data_file <- file.path(project_dir, "data", "simulations", subfolder_name, file_name)
  # df <- read.csv(data_file)
  # 
  # # Run all estimates on the DGP loaded into memory
  # output_df[startrow:endrow, ] <- simulation_estimates(
  #   i, df, sample, lin, hom, missing_type, missing_pct, conf
  # )
  
  ##########################################################
  # Third DGP
  dgp_num = 3; startrow = 1 + endrow; endrow = endrow + n_methods
  sample = 5000; lin = "lin"; hom = "hom"
  missing_type = "MCAR"; missing_pct = "98"; conf = "unconf"
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
  sample = 5000; lin = "lin"; hom = "hom"
  missing_type = "MCAR"; missing_pct = "98"; conf = "conf"
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
  sample = 5000; lin = "lin"; hom = "hom"
  missing_type = "MCAR"; missing_pct = "99"; conf = "unconf"
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
  sample = 5000; lin = "lin"; hom = "hom"
  missing_type = "MCAR"; missing_pct = "99"; conf = "conf"
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
}
close(output_file)

##########################################################
# Save simulation output

# Create new "snapshots" folder, if doesn't exist
snapshots_subfolder = file.path(project_dir, "outputs", "snapshots", datestamp)
ifelse(!dir.exists(snapshots_subfolder), dir.create(snapshots_subfolder), FALSE)

# Write overall simulation table for later access
save_simulation_output(output_df, "sim_", datestamp_input = datestamp, snapshot = TRUE)

# Collapse results into table of averages
mean_simulations <- summarize_simulation_results(output_df)

# Format table outputs
mean_simulations <- format_simulation_outputs(mean_simulations)

# Write overall simulation summary to a CSV file
save_simulation_output(mean_simulations, "summary_table_", 
                       datestamp_input = datestamp, snapshot = TRUE)

##########################################################
# Save as Latex tables for use in paper / presentation

### First output table
#  - Semi-supervised
#  - MCAR 90%
#  - IPW with logistic regression
#  - Linear Outcome
#  - Homogeneous Treatment Effect
xtable_df <- prepare_per_estimator_table(mean_simulations, "IPW")
# For presentation - with highlights
xtable_output <- format_per_estimator_xtable(xtable_df, n_sim, "MCAR", 
                                             highlight = TRUE)
save_per_estimator_to_latex(xtable_output, "IPW", datestamp_input = datestamp, 
                            snapshot = TRUE, highlight = TRUE)
save_per_estimator_to_latex(xtable_output, "IPW", datestamp_input = datestamp, 
                            snapshot = FALSE, highlight = TRUE)
# For paper - no highlights
xtable_output <- format_per_estimator_xtable(xtable_df, n_sim, "MCAR", 
                                             highlight = FALSE)
save_per_estimator_to_latex(xtable_output, "IPW", datestamp_input = datestamp, 
                            snapshot = TRUE, highlight = FALSE)
save_per_estimator_to_latex(xtable_output, "IPW", datestamp_input = datestamp, 
                            snapshot = FALSE, highlight = FALSE)

### Second output table
#  - Semi-supervised
#  - MCAR 90%
#  - IPW with BART
#  - Linear Outcome
#  - Homogeneous Treatment Effect
xtable_df <- prepare_per_estimator_table(mean_simulations, "IPW_BART")
# For presentation - with highlights
xtable_output <- format_per_estimator_xtable(xtable_df, n_sim, "MCAR", 
                                             highlight = TRUE)
save_per_estimator_to_latex(xtable_output, "IPW_BART", datestamp_input = datestamp, 
                            snapshot = TRUE, highlight = TRUE)
save_per_estimator_to_latex(xtable_output, "IPW_BART", datestamp_input = datestamp, 
                            snapshot = FALSE, highlight = TRUE)
# For paper - no highlights
xtable_output <- format_per_estimator_xtable(xtable_df, n_sim, "MCAR", 
                                             highlight = FALSE)
save_per_estimator_to_latex(xtable_output, "IPW_BART", datestamp_input = datestamp, 
                            snapshot = TRUE, highlight = FALSE)
save_per_estimator_to_latex(xtable_output, "IPW_BART", datestamp_input = datestamp, 
                            snapshot = FALSE, highlight = FALSE)

### Third output table
#  - Semi-supervised
#  - MCAR 90%
#  - IPW with BART
#  - Linear Outcome
#  - Homogeneous Treatment Effect
xtable_df <- prepare_per_estimator_table(mean_simulations, "TMLE")
# For presentation - with highlights
xtable_output <- format_per_estimator_xtable(xtable_df, n_sim, "MCAR", 
                                             highlight = TRUE)
save_per_estimator_to_latex(xtable_output, "TMLE", datestamp_input = datestamp, 
                            snapshot = TRUE, highlight = TRUE)
save_per_estimator_to_latex(xtable_output, "TMLE", datestamp_input = datestamp, 
                            snapshot = FALSE, highlight = TRUE)
# For paper - no highlights
xtable_output <- format_per_estimator_xtable(xtable_df, n_sim, "MCAR", 
                                             highlight = FALSE)
save_per_estimator_to_latex(xtable_output, "TMLE", datestamp_input = datestamp, 
                            snapshot = TRUE, highlight = FALSE)
save_per_estimator_to_latex(xtable_output, "TMLE", datestamp_input = datestamp, 
                            snapshot = FALSE, highlight = FALSE)

### Fourth output table
#  - Semi-supervised
#  - MCAR 90%
#  - BCF
#  - Linear Outcome
#  - Homogeneous Treatment Effect
xtable_df <- prepare_per_estimator_table(mean_simulations, "BCF")
# For presentation - with highlights
xtable_output <- format_per_estimator_xtable(xtable_df, n_sim, "MCAR", 
                                             highlight = TRUE)
save_per_estimator_to_latex(xtable_output, "BCF", datestamp_input = datestamp, 
                            snapshot = TRUE, highlight = TRUE)
save_per_estimator_to_latex(xtable_output, "BCF", datestamp_input = datestamp, 
                            snapshot = FALSE, highlight = TRUE)
# For paper - no highlights
xtable_output <- format_per_estimator_xtable(xtable_df, n_sim, "MCAR", 
                                             highlight = FALSE)
save_per_estimator_to_latex(xtable_output, "BCF", datestamp_input = datestamp, 
                            snapshot = TRUE, highlight = FALSE)
save_per_estimator_to_latex(xtable_output, "BCF", datestamp_input = datestamp, 
                            snapshot = FALSE, highlight = FALSE)

### Fifth output table
#  - Complete case
#  - MCAR 90%
#  - TMLE
#  - Linear Outcome
#  - Homogeneous Treatment Effect
xtable_df <- prepare_all_estimator_table(mean_simulations, "Complete Case")
# For presentation - with highlights
xtable_output <- format_all_estimator_xtable(xtable_df, n_sim, "MCAR", 
                                             highlight = TRUE)
save_all_estimator_to_latex(xtable_output, "CC", datestamp_input = datestamp, 
                            snapshot = TRUE, highlight = TRUE)
save_all_estimator_to_latex(xtable_output, "CC", datestamp_input = datestamp, 
                            snapshot = FALSE, highlight = TRUE)
# For paper - no highlights
xtable_output <- format_all_estimator_xtable(xtable_df, n_sim, "MCAR", 
                                             highlight = FALSE)
save_all_estimator_to_latex(xtable_output, "CC", datestamp_input = datestamp, 
                            snapshot = TRUE, highlight = FALSE)
save_all_estimator_to_latex(xtable_output, "CC", datestamp_input = datestamp, 
                            snapshot = FALSE, highlight = FALSE)

### Sixth output table
#  - Semi-supervised
#  - MCAR 90%
#  - TMLE
#  - Linear Outcome
#  - Homogeneous Treatment Effect
xtable_df <- prepare_all_estimator_table(mean_simulations, "Semi-supervised")
# For presentation - with highlights
xtable_output <- format_all_estimator_xtable(xtable_df, n_sim, "MCAR", 
                                             highlight = TRUE)
save_all_estimator_to_latex(xtable_output, "SSL", datestamp_input = datestamp, 
                            snapshot = TRUE, highlight = TRUE)
save_all_estimator_to_latex(xtable_output, "SSL", datestamp_input = datestamp, 
                            snapshot = FALSE, highlight = TRUE)
# For paper - no highlights
xtable_output <- format_all_estimator_xtable(xtable_df, n_sim, "MCAR", 
                                             highlight = FALSE)
save_all_estimator_to_latex(xtable_output, "SSL", datestamp_input = datestamp, 
                            snapshot = TRUE, highlight = FALSE)
save_all_estimator_to_latex(xtable_output, "SSL", datestamp_input = datestamp, 
                            snapshot = FALSE, highlight = FALSE)
