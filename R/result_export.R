# Preparing and formatting outputs for paper

##########################################################
# Save simulation output

# Write overall simulation table for later access
save_simulation_output <- function(simulation_results_df, prefix){
  csv_datetime <- paste0(prefix, format(Sys.time(), "%Y%m%d%H%M"), ".csv")
  output.file.name <- file.path(project_dir, "outputs", 
                                "simulation", csv_datetime)
  write.csv(simulation_results_df, output.file.name, row.names = FALSE)
}

####### Intended usage
# # Save the raw simulation results
# save_simulation_output(output_df, "sim_")
# # Save the collapsed simulation results
# save_simulation_output(mean_simulations, "summary_table_")

# Collapse results into table of averages
# Assumes column name structure
summarize_simulation_results <- function(simulation_results_df){
  mean_simulations <- simulation_results_df %>% 
    group_by(estimator, method, outcome, effect_type, 
             unlabeled_mechanism, pct_unlab, confoundedness, 
             other_notes, sample_size) %>% 
    summarise(true_tau = mean(true_tau), 
              ATE = mean(ATE), 
              ATE_CI_lb = mean(ATE_CI_lb), 
              ATE_CI_ub = mean(ATE_CI_ub), 
              rmse = mean(rmse), 
              cover = mean(cover)) %>%
    mutate(bias = ATE - true_tau)
  colnames(mean_simulations) <- c(
    "Estimator", "Method", "Outcome Model", "Treatment Effect Type", 
    "Unlabeled Mechanism", "% Unlabeled", "Confounding", 
    "Other_Notes", "N", "Tau", "ATE", "ATE_CI_lb", 
    "ATE_CI_ub", "RMSE", "Coverage", "Bias"
  )
  return(mean_simulations)
}
####### Intended usage
# # Collapse the simulation results
# mean_simulations <- summarize_simulation_results(output_df)

# Post-processing / formatting of table outputs
# Assumes column name structure
format_simulation_outputs <- function(sim_table_df){
  sim_table_df[, "Treatment"] = ifelse(
    sim_table_df$Confounding == "confounded", "Confounded", 
    "Experimental")
  sim_table_df[, "Propensities"] = ifelse(
    sim_table_df$Other_Notes == "True propensities", "True", 
    "Estimated")
  sim_table_df[, "Method"] = ifelse(
    sim_table_df$Method == "Complete Case", "CC", 
    "SSL")
  sim_table_df[, "# labeled"] = (1 - as.vector(sapply(sim_table_df[, '% Unlabeled'], as.numeric))/100)*sim_table_df$N
  sim_table_df$N = format(sim_table_df$N, digits = 0, big.mark = ",")
  sim_table_df[, "# labeled"] = format(as.vector(sapply(sim_table_df[, '# labeled'], as.numeric)), digits = 0, big.mark = ",")
  sim_table_df[, "Approach"] = ifelse(
    (sim_table_df$Other_Notes == "") & 
      (sim_table_df$Method == "CC"), "Complete Case", 
    ifelse((sim_table_df$Other_Notes == "") & 
             (sim_table_df$Method == "SSL"), "Semi-supervised", 
           ifelse((sim_table_df$Other_Notes == "True propensities") & 
                    (sim_table_df$Method == "CC"), "True propensities", "Other")))
  return(sim_table_df)
}
####### Intended usage
# # Format the collapsed simulation results
# mean_simulations <- format_simulation_outputs(mean_simulations)

##########################################################
## Save as Latex tables for use in paper / presentation

## First table functions: the "per-estimator" table
## (which reports 7 rows of results for each estimator)
prepare_per_estimator_table <- function(sim_table_df, estimator){
  xtable_df_1 <- sim_table_df[
    (sim_table_df$Estimator == estimator) & 
      (sim_table_df[, 'Outcome Model'] == "Linear") & 
      (sim_table_df[, 'Treatment'] == "Confounded") & 
      (sim_table_df[, 'Treatment Effect Type'] == "Homogeneous") & 
      (sim_table_df[, 'Unlabeled Mechanism'] == "MCAR") & 
      (sim_table_df[, 'Method'] == "CC") & 
      (sim_table_df[, '% Unlabeled'] == "90") & 
      (sim_table_df[, 'Propensities'] == "Estimated"), 
    c("Approach", "N", "# labeled", "ATE", "RMSE", "Bias", "Coverage")]
  xtable_df_2 <- sim_table_df[
    (sim_table_df$Estimator == estimator) & 
      (sim_table_df[, 'Outcome Model'] == "Linear") & 
      (sim_table_df[, 'Treatment'] == "Confounded") & 
      (sim_table_df[, 'Treatment Effect Type'] == "Homogeneous") & 
      (sim_table_df[, 'Unlabeled Mechanism'] == "MCAR") & 
      (sim_table_df[, 'Method'] == "SSL") & 
      (sim_table_df[, '% Unlabeled'] == "90") & 
      (sim_table_df[, 'Propensities'] == "Estimated"), 
    c("Approach", "N", "# labeled", "ATE", "RMSE", "Bias", "Coverage")]
  xtable_df_3 <- sim_table_df[
    (sim_table_df$Estimator == estimator) & 
      (sim_table_df[, 'Outcome Model'] == "Linear") & 
      (sim_table_df[, 'Treatment'] == "Confounded") & 
      (sim_table_df[, 'Treatment Effect Type'] == "Homogeneous") & 
      (sim_table_df[, 'Unlabeled Mechanism'] == "MCAR") & 
      (sim_table_df[, 'Method'] == "CC") & 
      (sim_table_df[, '% Unlabeled'] == "90") & 
      (sim_table_df[, 'Propensities'] == "True") & 
      (sim_table_df[, 'N'] == "5,000"), 
    c("Approach", "N", "# labeled", "ATE", "RMSE", "Bias", "Coverage")]
  xtable_df <- rbind(xtable_df_1, xtable_df_2, xtable_df_3)
  return(xtable_df)
}
####### Intended usage
# # Summarize simulation results into df ready for export
# xtable_df <- prepare_per_estimator_table(mean_simulations, "IPW")
# xtable_df <- prepare_per_estimator_table(mean_simulations, "IPW_BART")
# xtable_df <- prepare_per_estimator_table(mean_simulations, "TMLE")
# xtable_df <- prepare_per_estimator_table(mean_simulations, "BCF")

format_per_estimator_xtable <- function(xtable_df, n_sim, 
                                        missing_type, highlight = TRUE){
  if (highlight) {
    alignment <- c("l", "l", "l", "l", "c", "c", 
                   ">{\\columncolor[gray]{.9}}c", ">{\\columncolor[gray]{.9}}c")
  } else {
    alignment <- c("l", "l", "l", "l", "c", "c", "c", "c")
  }
  
  output <- xtable(
    xtable_df, type = "latex", 
    align = alignment, 
    caption = paste0(n_sim, " simulations with outcomes ", missing_type)
  )
  return(output)
}
####### Intended usage
# # Format simulation result table as xtable output
# xtable_output <- format_per_estimator_xtable(xtable_df, n_sim, "MCAR")

save_per_estimator_to_latex <- function(xtable_output, estimator, 
                                        highlight = TRUE){
  if (highlight){
    file_name = paste0("outputs/simulations_", estimator, 
                       "_LIN_HOM_MCAR_90_HIGHLIGHT.tex")
    print(xtable_output,
          floating = TRUE, latex.environments = "center",
          include.rownames = FALSE, size="small", 
          add.to.row=list(list(3,4,5,6),c("\\rowcolor[gray]{.8} ", 
                                          "\\rowcolor[gray]{.7} ", 
                                          "\\rowcolor[gray]{.6} ", 
                                          "\\rowcolor[gray]{.5} ")), 
          file = file_name
    )
  } else {
    file_name = paste0("outputs/simulations_", estimator, 
                       "_LIN_HOM_MCAR_90.tex")
    print(xtable_output,
          floating = TRUE, latex.environments = "center",
          include.rownames = FALSE, size="small", 
          file = file_name
    )
  }
}
####### Intended usage
# # Save result as .tex file
# save_per_estimator_to_latex(xtable_output, "IPW")
# save_per_estimator_to_latex(xtable_output, "IPW_BART")
# save_per_estimator_to_latex(xtable_output, "TMLE")
# save_per_estimator_to_latex(xtable_output, "BCF")

## First table functions: the "all estimator" table
## (which reports 12 rows of results for a given approach: ssl, cc, etc..)
prepare_all_estimator_table <- function(sim_table_df, approach){
  xtable_df_1 <- sim_table_df[
    (sim_table_df$Estimator == "IPW") & 
      (sim_table_df[, 'Outcome Model'] == "Linear") & 
      (sim_table_df[, 'Treatment'] == "Experimental") & 
      (sim_table_df[, 'Treatment Effect Type'] == "Homogeneous") & 
      (sim_table_df[, 'Unlabeled Mechanism'] == "MCAR") & 
      (sim_table_df[, 'Approach'] == approach) & 
      (sim_table_df[, '% Unlabeled'] == "90"), 
    c("Estimator", "N", "# labeled", "ATE", "RMSE", "Bias", "Coverage")]
  xtable_df_2 <- sim_table_df[
    (sim_table_df$Estimator == "IPW_BART") & 
      (sim_table_df[, 'Outcome Model'] == "Linear") & 
      (sim_table_df[, 'Treatment'] == "Experimental") & 
      (sim_table_df[, 'Treatment Effect Type'] == "Homogeneous") & 
      (sim_table_df[, 'Unlabeled Mechanism'] == "MCAR") & 
      (sim_table_df[, 'Approach'] == approach) & 
      (sim_table_df[, '% Unlabeled'] == "90"), 
    c("Estimator", "N", "# labeled", "ATE", "RMSE", "Bias", "Coverage")]
  xtable_df_3 <- sim_table_df[
    (sim_table_df$Estimator == "TMLE") & 
      (sim_table_df[, 'Outcome Model'] == "Linear") & 
      (sim_table_df[, 'Treatment'] == "Experimental") & 
      (sim_table_df[, 'Treatment Effect Type'] == "Homogeneous") & 
      (sim_table_df[, 'Unlabeled Mechanism'] == "MCAR") & 
      (sim_table_df[, 'Approach'] == approach) & 
      (sim_table_df[, '% Unlabeled'] == "90"), 
    c("Estimator", "N", "# labeled", "ATE", "RMSE", "Bias", "Coverage")]
  xtable_df_4 <- sim_table_df[
    (sim_table_df$Estimator == "BCF") & 
      (sim_table_df[, 'Outcome Model'] == "Linear") & 
      (sim_table_df[, 'Treatment'] == "Experimental") & 
      (sim_table_df[, 'Treatment Effect Type'] == "Homogeneous") & 
      (sim_table_df[, 'Unlabeled Mechanism'] == "MCAR") & 
      (sim_table_df[, 'Approach'] == approach) & 
      (sim_table_df[, '% Unlabeled'] == "90"), 
    c("Estimator", "N", "# labeled", "ATE", "RMSE", "Bias", "Coverage")]
  xtable_df <- rbind(xtable_df_1, xtable_df_2, 
                     xtable_df_3, xtable_df_4)
  Estimator_Old <- as.character(xtable_df$Estimator)
  xtable_df$Estimator <- ifelse(
    (Estimator_Old == "IPW_BART"), "IPW (BART)",
    ifelse((Estimator_Old == "IPW"), "IPW (logistic)",
           Estimator_Old)
  )
  return(xtable_df)
}
####### Intended usage
# # Summarize simulation results into df ready for export
# xtable_df <- prepare_all_estimator_table(mean_simulations, "Semi-supervised")

format_all_estimator_xtable <- function(xtable_df, n_sim, missing_type, 
                                        highlight = TRUE){
  if (highlight) {
    alignment <- c("l", "l", "l", "l", "c", "c", 
                   ">{\\columncolor[gray]{.9}}c", 
                   ">{\\columncolor[gray]{.9}}c")
  } else {
    alignment <- c("l", "l", "l", "l", "c", "c", "c", "c")
  }
  
  output <- xtable(
    xtable_df, type = "latex", 
    align = alignment, 
    caption = paste0(n_sim, " simulations with outcomes ", missing_type)
  )
  return(output)
}
####### Intended usage
# # Format simulation result table as xtable output
# xtable_output <- format_all_estimator_xtable(xtable_df, n_sim, "MCAR")

save_all_estimator_to_latex <- function(xtable_output, approach, 
                                        highlight = TRUE){
  if (highlight) {
    file_name = paste0("outputs/simulations_", approach, 
                       "_ALL_EXP_LIN_HOM_MCAR_90_HIGHLIGHT.tex")
  }
  else {
    file_name = paste0("outputs/simulations_", approach, 
                       "_ALL_EXP_LIN_HOM_MCAR_90.tex")
  }
  print(xtable_output,
        floating = TRUE, latex.environments = "center",
        include.rownames = FALSE, size="small", 
        file = file_name
  )
}
####### Intended usage
# # Save result as .tex file
# save_all_estimator_to_latex(xtable_output, "SSL")
