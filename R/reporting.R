## Storing / reporting program results

create_simulation_row <- function(
  i, estimator, method, outcome, tau, unlabeled_mech, 
  pct_unlab, confoundedness, other_notes, sample_size, 
  simulation_object
){
  return(c(list(i, estimator, method, outcome, tau, unlabeled_mech, 
              pct_unlab, confoundedness, other_notes, sample_size), 
              as.list(simulation_object)))
}

model_run <- function(
  sim_num, df, sample = 5000, lin = "lin", hom = "hom", 
  missing_type = "MCAR", missing_pct = "90", conf = "unconf", 
  type="ssl", alpha = 0.05, method = "IPW", other_notes, 
  estimate_propensity = TRUE
){
  # Generate estimate
  if (method == "IPW"){
    estimate = IPW_linear(df, type=type, alpha = alpha, 
                          estimate_propensity = estimate_propensity)
  } else if (method == "IPW_BART"){
      estimate = IPW_BART(df, type=type, alpha = alpha, 
                            estimate_propensity = estimate_propensity)
  } else if (method == "IPW_XBART"){
    estimate = IPW_XBART(df, type=type, alpha = alpha, 
                        estimate_propensity = estimate_propensity)
  } else if (method == "TMLE"){
    estimate = TMLE_estimator(df, type=type, 
                              estimate_propensity = estimate_propensity)
  } else if (method == "TMLE_XBART"){
    estimate = TMLE_XBART_estimator(df, type=type, 
                              estimate_propensity = estimate_propensity)
  } else if (method == "BCF"){
    estimate = BCF_estimator(df, type=type, alpha = alpha, 
                             estimate_propensity = estimate_propensity)
  } else if (method == "BCF_XBART"){
    estimate = BCF_XBART_estimator(df, type=type, alpha = alpha, 
                             estimate_propensity = estimate_propensity)
  } else if (method == "XBCF"){
    estimate = XBCF_estimator(df, type=type, alpha = alpha, 
                             estimate_propensity = estimate_propensity)
  }
  
  # Format outputs for table
  estimate_type = ifelse(type == "ssl", "SSL", ifelse(type == "complete_case", "Complete Case", "Full Data"))
  linear_nonlinear = ifelse(lin == "lin", "Linear", "Nonlinear")
  effect_type = ifelse(hom == "hom", "Homogeneous", "Heterogeneous")
  confoundedness = ifelse(conf == "conf", "confounded", "unconfounded")
  
  return(create_simulation_row(
    sim_num, method, estimate_type, linear_nonlinear, effect_type, 
    missing_type, missing_pct, confoundedness, other_notes, sample, 
    estimate
  ))
}

simulation_estimates <- function(sim_num, df, 
  # Details about the DGP
  sample = 5000, lin = "lin", hom = "hom", missing_type = "MCAR", 
  missing_pct = "90", conf = "unconf"
){
  # Allocate space for outputs
  row_num = 1
  num_rows = 12
  simulation_result_output = data.frame(
    sim_num=rep(NA, num_rows), 
    estimator=rep("", num_rows), 
    method=rep("", num_rows), 
    outcome=rep("", num_rows), 
    effect_type=rep("", num_rows), 
    unlabeled_mechanism=rep("", num_rows), 
    pct_unlab=rep("", num_rows), 
    confoundedness=rep("", num_rows), 
    other_notes=rep("", num_rows), 
    sample_size=rep(NA, num_rows), 
    true_tau=rep(NA, num_rows), 
    ATE=rep(NA, num_rows), 
    ATE_CI_lb=rep(NA, num_rows), 
    ATE_CI_ub=rep(NA, num_rows), 
    rmse=rep(NA, num_rows), 
    cover=rep(NA, num_rows), 
    stringsAsFactors=FALSE
  )
  
  ###############################################
  # IPW-Linear
  # Simple IPW estimator
  # Semi-supervised case
  simulation_result_output[row_num, ] <- model_run(
    sim_num, df, sample, lin, hom, missing_type, 
    missing_pct, conf, type="ssl", alpha = 0.05, method = "IPW", 
    other_notes="", estimate_propensity = TRUE
  )
  
  # Simple IPW estimator
  # Complete case analysis
  row_num = row_num + 1
  simulation_result_output[row_num, ] <- model_run(
    sim_num, df, sample, lin, hom, missing_type, 
    missing_pct, conf, type="complete_case", alpha = 0.05, method = "IPW", 
    other_notes="", estimate_propensity = TRUE
  )
  
  # Simple IPW estimator
  # Semi-supervised case
  # row_num = row_num + 1
  # simulation_result_output[row_num, ] <- model_run(
  #   sim_num, df, sample, lin, hom, missing_type, 
  #   missing_pct, conf, type="ssl", alpha = 0.05, method = "IPW", 
  #   other_notes="True propensities", estimate_propensity = FALSE
  # )
  
  # Simple IPW estimator
  # Complete case analysis
  row_num = row_num + 1
  simulation_result_output[row_num, ] <- model_run(
    sim_num, df, sample, lin, hom, missing_type, 
    missing_pct, conf, type="complete_case", alpha = 0.05, method = "IPW", 
    other_notes="True propensities", estimate_propensity = FALSE
  )
  
  ###############################################
  # TMLE
  # Run TMLE estimator with BART for prognostic and 
  # propensity models
  # Semi-supervised case
  row_num = row_num + 1
  simulation_result_output[row_num, ] <- model_run(
    sim_num, df, sample, lin, hom, missing_type, 
    missing_pct, conf, type="ssl", alpha = 0.05, method = "TMLE", 
    other_notes="", estimate_propensity = TRUE
  )

  # Run TMLE estimator with BART for prognostic and 
  # propensity models
  # Complete case analysis
  row_num = row_num + 1
  simulation_result_output[row_num, ] <- model_run(
    sim_num, df, sample, lin, hom, missing_type, 
    missing_pct, conf, type="complete_case", alpha = 0.05, method = "TMLE", 
    other_notes="", estimate_propensity = TRUE
  )

  # Run TMLE estimator with BART for prognostic and 
  # propensity models
  # Semi-supervised case
  # row_num = row_num + 1
  # simulation_result_output[row_num, ] <- model_run(
  #   sim_num, df, sample, lin, hom, missing_type, 
  #   missing_pct, conf, type="ssl", alpha = 0.05, method = "TMLE", 
  #   other_notes="True propensities", estimate_propensity = FALSE
  # )

  # Run TMLE estimator with BART for prognostic and 
  # propensity models
  # Complete case analysis
  row_num = row_num + 1
  simulation_result_output[row_num, ] <- model_run(
    sim_num, df, sample, lin, hom, missing_type, 
    missing_pct, conf, type="complete_case", alpha = 0.05, method = "TMLE", 
    other_notes="True propensities", estimate_propensity = FALSE
  )
  
  ###############################################
  # BCF
  # BCF estimator
  # Semi-supervised case
  row_num = row_num + 1
  simulation_result_output[row_num, ] <- model_run(
    sim_num, df, sample, lin, hom, missing_type, 
    missing_pct, conf, type="ssl", alpha = 0.05, method = "BCF", 
    other_notes="", estimate_propensity = TRUE
  )
  
  # BCF estimator
  # Complete case
  row_num = row_num + 1
  simulation_result_output[row_num, ] <- model_run(
    sim_num, df, sample, lin, hom, missing_type, 
    missing_pct, conf, type="complete_case", alpha = 0.05, method = "BCF", 
    other_notes="", estimate_propensity = TRUE
  )

  # BCF estimator
  # Semi-supervised
  # True propensities
  # row_num = row_num + 1
  # simulation_result_output[row_num, ] <- model_run(
  #   sim_num, df, sample, lin, hom, missing_type, 
  #   missing_pct, conf, type="ssl", alpha = 0.05, method = "BCF", 
  #   other_notes="True propensities", estimate_propensity = FALSE
  # )
  
  # BCF estimator
  # Complete case
  # True propensities
  row_num = row_num + 1
  simulation_result_output[row_num, ] <- model_run(
    sim_num, df, sample, lin, hom, missing_type, 
    missing_pct, conf, type="complete_case", alpha = 0.05, method = "BCF", 
    other_notes="True propensities", estimate_propensity = FALSE
  )
  
  ###############################################
  # IPW-BART
  # IPW-BART estimator
  # Semi-supervised case
  row_num = row_num + 1
  simulation_result_output[row_num, ] <- model_run(
    sim_num, df, sample, lin, hom, missing_type, 
    missing_pct, conf, type="ssl", alpha = 0.05, method = "IPW_BART", 
    other_notes="", estimate_propensity = TRUE
  )
  
  # IPW-BART estimator
  # Complete case analysis
  row_num = row_num + 1
  simulation_result_output[row_num, ] <- model_run(
    sim_num, df, sample, lin, hom, missing_type, 
    missing_pct, conf, type="complete_case", alpha = 0.05, method = "IPW_BART", 
    other_notes="", estimate_propensity = TRUE
  )
  
  # IPW-BART estimator
  # Semi-supervised case
  # row_num = row_num + 1
  # simulation_result_output[row_num, ] <- model_run(
  #   sim_num, df, sample, lin, hom, missing_type, 
  #   missing_pct, conf, type="ssl", alpha = 0.05, method = "IPW_BART", 
  #   other_notes="True propensities", estimate_propensity = FALSE
  # )
  
  # IPW-BART estimator
  # Complete case analysis
  row_num = row_num + 1
  simulation_result_output[row_num, ] <- model_run(
    sim_num, df, sample, lin, hom, missing_type, 
    missing_pct, conf, type="complete_case", alpha = 0.05, method = "IPW_BART", 
    other_notes="True propensities", estimate_propensity = FALSE
  )

  return(simulation_result_output)
}

unpack_simulation_results <- function(
  i, row_num, estimator, method, outcome, tau, unlabeled_mech, 
  pct_unlab, confoundedness, other_notes, sample_size, 
  output_object, simulation_object
){
  output_object[row_num, 1] = i; output_object[row_num, 2] = estimator
  output_object[row_num, 3] = method; output_object[row_num, 4] = outcome
  output_object[row_num, 5] = tau; output_object[row_num, 6] = unlabeled_mech
  output_object[row_num, 7] = pct_unlab; output_object[row_num, 8] = confoundedness
  output_object[row_num, 9] = other_notes
  output_object[row_num, 10] = sample_size
  output_object[row_num, 11:16] = simulation_object
  return(output_object)
}


propensity.pairplot <- function(tau.hat, my_data){
  # Generate and save pairplot of predicted propensities vs 
  # actual using variety of propensity modeling techniques
  # @param tau.hat: list of estimated treatment probabilities
  # @param my_data: list of components of a dataset
  # @return: pairplot of propensity scores
  curr_datetime <- format(Sys.time(), format="%Y%m%d_%H%M%S")
  pairplot_filename <- paste0("outputs/pairplot ",curr_datetime,".png")
  png(pairplot_filename)
  pairs(cbind(Correct.Logit = tau.hat$W.hat.logit.correct, 
              General.Logit = tau.hat$W.hat.logit.kitchen.sink, 
              LASSO = tau.hat$W.hat.lasso, 
              Ridge = tau.hat$W.hat.ridge, 
              GRF = tau.hat$W.hat.grf, 
              XBART = tau.hat$W.hat.xbart, 
              Tau = my_data$T.prob))
  dev.off()
}