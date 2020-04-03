# ========================================================
# Runner script for parallel simulations
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

# Create new "snapshots" folder, if doesn't exist
snapshots_subfolder = file.path(project_dir, "outputs", "snapshots", datestamp)
ifelse(!dir.exists(snapshots_subfolder), dir.create(snapshots_subfolder), FALSE)

# Simulation DGPs
DGP.list <- data.frame(
  sample = c(5000, 5000, 5000, 5000, 5000, 5000),
  lin = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
  hom = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
  missing_type = c("MCAR", "MCAR", "MCAR", "MCAR", "MCAR", "MCAR"),
  missing_pct = c(0.90, 0.90, 0.98, 0.98, 0.99, 0.99),
  conf = c(FALSE, TRUE, FALSE, TRUE, FALSE, TRUE),
  stringsAsFactors = FALSE
)

# Create outputs file to store results
n_methods = 12
n_sim = 1
n_dgps = nrow(DGP.list)
n_df = n_sim*n_dgps*n_methods
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

# Register parallel backend
numCores <- detectCores()
registerDoParallel(numCores)

# Run simulations, with the DGPs run in parallel
for (i in 1:n_sim){
  t = proc.time()[3]
  
  # Rows to save to `output_df`
  startrow = 1 + (i-1)*n_dgps*n_methods
  endrow = startrow + n_dgps*n_methods - 1
  
  # Run each of the simulated DGPs in parallel
  system.time(sim_iter_result <- foreach (dgp = 1:n_dgps, .combine=rbind) %dopar% {
    # Extract simulation directions
    n <- DGP.list[dgp, 'sample']
    linear <- DGP.list[dgp, 'lin']
    homogeneous <- DGP.list[dgp, 'hom']
    missing_mechanism <- DGP.list[dgp, 'missing_type']
    missing_pct <- DGP.list[dgp, 'missing_pct']
    confounded <- DGP.list[dgp, 'conf']
    
    # Set seed for each repetition of simulation
    random_seed_incr <- as.integer(paste0(dgp, paste0(rep("0", 4-length(paste0(dgp,i))), collapse=""), i))
    set.seed(random_seed_incr)
    
    # Generate BCF data
    dummy_data <- bcf_simulations(n = n, linear = linear, 
                                  homogeneous = homogeneous,
                                  missing_mechanism = missing_mechanism, 
                                  missing_pct = missing_pct, 
                                  confounded = confounded)
    causal_data = data.frame(cbind(Y = dummy_data$Y, 
                                   Z = dummy_data$Z, 
                                   dummy_data$X, 
                                   M = dummy_data$M, 
                                   tau = dummy_data$tau, 
                                   pi_x = dummy_data$pi_x))
    
    # Convert simulation instructions into formatted
    # directions for estimation and compilation of results
    lin = ifelse(linear, "lin", "nonlin")
    hom = ifelse(homogeneous, "hom", "het")
    missing = paste(floor(missing_pct*100))
    conf = ifelse(confounded, "conf", "unconf")

    # Run each of the `n_methods` estimators on `causal_data`
    # simulated dataset
    simulation_estimates(
      i, causal_data, n, lin, hom, missing_mechanism, missing, conf
    )
  })
  
  # Save each simulation iteration to `output_df`
  output_df[startrow:endrow, ] <- sim_iter_result
  
  # Write inidividual simulation results for intermediate access
  save_simulation_output(sim_iter_result, paste0("intermediate_sim_", i, "_"), datestamp_input = datestamp, snapshot = TRUE)
  
  t = proc.time()[3] - t
  cat("==================================================\n", file = output_file)
  cat(paste("Run time:", t, "\n"), file = output_file)
  cat(paste("Sim", i, "of", n_sim, "complete\n\n"), file = output_file)
}
close(output_file)

##########################################################
# Save simulation output

# Write overall simulation table for later access
save_simulation_output(output_df, "simulations_", datestamp_input = datestamp, snapshot = TRUE)

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
