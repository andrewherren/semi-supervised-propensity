sample_split <- function(K, N){
  # Return split indices based on N and K
  return(sample(1:K, N, replace = T))
}

tau_fold <- function(Y, Z, Z.est, Y.est){
  prop.error = Z - Z.est
  prog_error = Y - Y.est
  n.fold = length(Z)
  tau.est.full = ((sum(prop.error*Z)/n.fold)^(-1))*(sum(prop.error*prog_error)/n.fold)
  return(tau.est.full)
}

outcome_estimates <- function(X, Y, indices, fold){
  prognostic.model = bart(X[indices != fold, ], Y[indices != fold], keeptrees = T)
  return(colMeans(predict(prognostic.model, X[indices == fold, ])))
}

propensity_estimates <- function(X, Z, indices, fold){
  propensity.model = bart(X[indices != fold, ], Z[indices != fold], keeptrees = T)
  return(colMeans(pnorm(predict(propensity.model, X[indices == fold, ]))))
}

IPW_linear <- function(df, type="full", alpha = 0.05, 
                       estimate_propensity = TRUE){
  # Requires that the dataframe have the following columns:
  # Y, Z, x_1, x_2, x_3, x_4, x_5
  outcome_col <- "Y"
  treatment_col <- "Z"
  covariates <- c("x_1", "x_2", "x_3", "x_4", "x_5")
  
  if (type == "full"){
    propensity.data <- df
    outcome.data <- df
    N = length(df$Y)
  } else if (type == "ssl"){
    propensity.data <- df
    outcome.data <- df[df$M == 0, ]
    N = length(df[df$M == 0, "Y"])
  } else if (type == "complete_case"){
    propensity.data <- df[df$M == 0, ]
    outcome.data <- df[df$M == 0, ]
    N = length(df[df$M == 0, "Y"])
  }
  
  # Estimate propensities
  if (estimate_propensity){
    propensity.model <- glm(Z ~ x_1 + x_2 + x_3 + x_4 + x_5, 
                            family = binomial(link = "logit"), 
                            data = propensity.data)
    propensity.scores <- predict(propensity.model, 
                                 newdata = outcome.data[, covariates], 
                                 type = "response")
  } else{
    propensity.scores = outcome.data$pi_x
  }

  # Estimate ATE
  ATE.summand = ((outcome.data$Z - propensity.scores)*outcome.data$Y)/(propensity.scores*(1 - propensity.scores))
  ATE = sum(ATE.summand)/N
  
  # Estimate variance of ATE
  score = (outcome.data$Z - propensity.scores)*outcome.data[, covariates]
  ATE.df = data.frame(ATE.summand = ATE.summand, score = score)
  ATE.score.model <- lm(ATE.summand ~ ., data=ATE.df)
  SE = sqrt(sum((ATE.score.model$residuals)^2)/N)/sqrt(N)
  
  # Return estimand and confidence interval
  normal.cutoff <- qnorm(1-(alpha/2))
  ATE_CI_lb = ATE - normal.cutoff*SE
  ATE_CI_ub = ATE + normal.cutoff*SE
  
  # Estimate RMSE and CI coverage
  true_tau <- mean(outcome.data[, "tau"])
  rmse <- sqrt(mean((ATE - outcome.data[, "tau"])^2))
  cover <- (true_tau >= ATE_CI_lb) & (true_tau <= ATE_CI_ub)
  return(c(true_tau, ATE, ATE_CI_lb, ATE_CI_ub, rmse, cover))
}

TMLE_estimator <- function(df, type = "full", 
                           estimate_propensity = TRUE){
  # Requires that the dataframe have the following columns:
  # Y, Z, x_1, x_2, x_3, x_4, x_5
  outcome_col <- "Y"
  treatment_col <- "Z"
  covariates <- c("x_1", "x_2", "x_3", "x_4", "x_5")
  
  if (type == "full"){
    propensity.data <- df
    outcome.data <- df
    N = length(df$Y)
  } else if (type == "ssl"){
    propensity.data <- df
    outcome.data <- df[df$M == 0, ]
    N = length(df[df$M == 0, "Y"])
  } else if (type == "complete_case"){
    propensity.data <- df[df$M == 0, ]
    outcome.data <- df[df$M == 0, ]
    N = length(df[df$M == 0, "Y"])
  }
  
  W_a = outcome.data[, c(treatment_col, covariates)]
  W_1 = cbind(Z=1, outcome.data[, covariates])
  W_0 = cbind(Z=0, outcome.data[, covariates])
  Y = outcome.data[, outcome_col]

  # Fit prognostic model
  prognostic.model = bart(W_a, Y, keeptrees = T, verbose = F)
  # Q_a_pred = colMeans(predict(prognostic.model, W_a))
  Q_1_pred = colMeans(predict(prognostic.model, W_1))
  Q_0_pred = colMeans(predict(prognostic.model, W_0))
  Q_pred = cbind(Q_0_pred, Q_1_pred)
  
  # Build propensity model and predict propensity scores
  X.train = propensity.data[, covariates]
  Z.train = propensity.data[, treatment_col]
  X.pred = outcome.data[, covariates]
  Z.pred = outcome.data[, treatment_col]
  
  # Estimate tau and 95% confidence interval
  if (estimate_propensity){
    propensity.model = bart(X.train, Z.train, keeptrees = T, verbose = F)
    propensity.pred = colMeans(pnorm(predict(propensity.model, X.pred)))
  } else{
    propensity.pred = outcome.data$pi_x
  }
  
  tau.model <- tmle(Y, Z.pred, X.pred, 
                    Q = Q_pred, g1W = propensity.pred)
  ATE = tau.model$estimates$ATE$psi
  ATE_CI_lb = tau.model$estimates$ATE$CI[1]
  ATE_CI_ub = tau.model$estimates$ATE$CI[2]
  
  # Estimate RMSE and CI coverage
  true_tau <- mean(outcome.data[, "tau"])
  rmse <- sqrt(mean((ATE - outcome.data[, "tau"])^2))
  cover <- (true_tau >= ATE_CI_lb) & (true_tau <= ATE_CI_ub)
  return(c(true_tau, ATE, ATE_CI_lb, ATE_CI_ub, rmse, cover))
}

BCF_estimator <- function(df, type = "full", alpha = 0.05, 
                          estimate_propensity = TRUE){
  # Requires that the dataframe have the following columns:
  # Y, Z, x_1, x_2, x_3, x_4, x_5
  outcome_col <- "Y"
  treatment_col <- "Z"
  covariates <- c("x_1", "x_2", "x_3", "x_4", "x_5")
  
  if (type == "full"){
    propensity.data <- df
    outcome.data <- df
    N = length(df$Y)
  } else if (type == "ssl"){
    propensity.data <- df
    outcome.data <- df[df$M == 0, ]
    N = length(df[df$M == 0, "Y"])
  } else if (type == "complete_case"){
    propensity.data <- df[df$M == 0, ]
    outcome.data <- df[df$M == 0, ]
    N = length(df[df$M == 0, "Y"])
  }
  
  # Build propensity model and predict propensity scores
  X.train = propensity.data[, covariates]
  Z.train = propensity.data[, treatment_col]
  X.pred = outcome.data[, covariates]
  Z.pred = outcome.data[, treatment_col]
  
  # Estimate tau and 95% confidence interval
  if (estimate_propensity){
    propensity.model = bart(X.train, Z.train, keeptrees = T, verbose = F)
    propensity.pred = colMeans(pnorm(predict(propensity.model, X.pred)))
  } else{
    propensity.pred = outcome.data$pi_x
  }
  
  bcf.model <- bcf(outcome.data$Y, outcome.data$Z, as.matrix(X.pred), 
                   pihat = propensity.pred, nburn=100, nsim=100)
  ATE = mean(colMeans(bcf.model$tau))
  ATE_CI_lb = unname(quantile(colMeans(bcf.model$tau), (alpha/2)))
  ATE_CI_ub = unname(quantile(colMeans(bcf.model$tau), (1-(alpha/2))))
  
  # Estimate RMSE and CI coverage
  true_tau <- mean(outcome.data[, "tau"])
  rmse <- sqrt(mean((colMeans(bcf.model$tau) - outcome.data[, "tau"])^2))
  cover <- (true_tau >= ATE_CI_lb) & (true_tau <= ATE_CI_ub)
  return(c(true_tau, ATE, ATE_CI_lb, ATE_CI_ub, rmse, cover))
}

IPW_BART <- function(df, type="full", alpha = 0.05, 
                       estimate_propensity = TRUE){
  # Requires that the dataframe have the following columns:
  # Y, Z, x_1, x_2, x_3, x_4, x_5
  outcome_col <- "Y"
  treatment_col <- "Z"
  covariates <- c("x_1", "x_2", "x_3", "x_4", "x_5")
  
  if (type == "full"){
    propensity.data <- df
    outcome.data <- df
    N = length(df$Y)
  } else if (type == "ssl"){
    propensity.data <- df
    outcome.data <- df[df$M == 0, ]
    N = length(df[df$M == 0, "Y"])
  } else if (type == "complete_case"){
    propensity.data <- df[df$M == 0, ]
    outcome.data <- df[df$M == 0, ]
    N = length(df[df$M == 0, "Y"])
  }
  
  # Build propensity model and predict propensity scores
  X.train = propensity.data[, covariates]
  Z.train = propensity.data[, treatment_col]
  X.pred = outcome.data[, covariates]
  Z.pred = outcome.data[, treatment_col]
  
  # Estimate propensities
  if (estimate_propensity){
    propensity.model = bart(X.train, Z.train, keeptrees = T, verbose = F)
    propensity.scores = colMeans(pnorm(predict(propensity.model, X.pred)))
  } else{
    propensity.scores = outcome.data$pi_x
  }
  
  # Estimate ATE
  ATE.summand = ((outcome.data$Z - propensity.scores)*outcome.data$Y)/(propensity.scores*(1 - propensity.scores))
  ATE = sum(ATE.summand)/N
  
  # Estimate variance of ATE
  score = (outcome.data$Z - propensity.scores)*outcome.data[, covariates]
  ATE.df = data.frame(ATE.summand = ATE.summand, score = score)
  ATE.score.model <- lm(ATE.summand ~ ., data=ATE.df)
  SE = sqrt(sum((ATE.score.model$residuals)^2)/N)/sqrt(N)
  
  # Return estimand and confidence interval
  normal.cutoff <- qnorm(1-(alpha/2))
  ATE_CI_lb = ATE - normal.cutoff*SE
  ATE_CI_ub = ATE + normal.cutoff*SE
  
  # Estimate RMSE and CI coverage
  true_tau <- mean(outcome.data[, "tau"])
  rmse <- sqrt(mean((ATE - outcome.data[, "tau"])^2))
  cover <- (true_tau >= ATE_CI_lb) & (true_tau <= ATE_CI_ub)
  return(c(true_tau, ATE, ATE_CI_lb, ATE_CI_ub, rmse, cover))
}

IPW_XBART <- function(df, type="full", alpha = 0.05, 
                     estimate_propensity = TRUE){
  # Requires that the dataframe have the following columns:
  # Y, Z, x_1, x_2, x_3, x_4, x_5
  outcome_col <- "Y"
  treatment_col <- "Z"
  covariates <- c("x_1", "x_2", "x_3", "x_4", "x_5")
  
  if (type == "full"){
    propensity.data <- df
    outcome.data <- df
    N = length(df$Y)
  } else if (type == "ssl"){
    propensity.data <- df
    outcome.data <- df[df$M == 0, ]
    N = length(df[df$M == 0, "Y"])
  } else if (type == "complete_case"){
    propensity.data <- df[df$M == 0, ]
    outcome.data <- df[df$M == 0, ]
    N = length(df[df$M == 0, "Y"])
  }
  
  # Build propensity model and predict propensity scores
  X.train = propensity.data[, covariates]
  Z.train = propensity.data[, treatment_col]
  X.pred = outcome.data[, covariates]
  Z.pred = outcome.data[, treatment_col]
  
  # Estimate propensities
  if (estimate_propensity){
    
    get_XBART_params <- function(y, p){
      XBART_params = list(M = 15,
                          L = 1,
                          nsweeps = 150,
                          Nmin = 100,
                          alpha = 0.95,
                          beta = 1.25,
                          mtry = p,
                          burnin = 15)
      num_trees = XBART_params$M
      XBART_params$max_depth = 250
      XBART_params$Ncutpoints = 100;
      XBART_params$tau = var(y)/(num_trees)
      XBART_params$a = 0.000001; XBART_params$b = 0.000001;
      return(XBART_params)
    }
    params = get_XBART_params(Z.train, ncol(X.train))
    dcat = 0
    parl = F
    
    propensity.model = XBART.Probit(
      as.matrix(Z.train), as.matrix(X.train), as.matrix(X.train), num_trees = params$M,
      L = 1, num_sweeps = params$nsweeps, 
      max_depth = params$max_depth, Nmin = 10,
      num_cutpoints = params$Ncutpoints,
      alpha = params$alpha, beta = params$beta,
      tau = params$tau, s= 1,kap = 1,
      mtry = params$mtry, p_categorical = dcat,
      draw_sigma = FALSE, m_update_sigma = TRUE,
      parallel = parl, random_seed = 10
    )
    propensity.scores = apply(
      pnorm(predict(propensity.model, as.matrix(X.pred))[,params$burnin:params$nsweeps]), 1, mean
    )
    # 
    # propensity.model = bart(X.train, Z.train, keeptrees = T, verbose = F)
    # propensity.scores = colMeans(pnorm(predict(propensity.model, X.pred)))
  } else{
    propensity.scores = outcome.data$pi_x
  }
  
  # Estimate ATE
  ATE.summand = ((outcome.data$Z - propensity.scores)*outcome.data$Y)/(propensity.scores*(1 - propensity.scores))
  ATE = sum(ATE.summand)/N
  
  # Estimate variance of ATE
  score = (outcome.data$Z - propensity.scores)*outcome.data[, covariates]
  ATE.df = data.frame(ATE.summand = ATE.summand, score = score)
  ATE.score.model <- lm(ATE.summand ~ ., data=ATE.df)
  SE = sqrt(sum((ATE.score.model$residuals)^2)/N)/sqrt(N)
  
  # Return estimand and confidence interval
  normal.cutoff <- qnorm(1-(alpha/2))
  ATE_CI_lb = ATE - normal.cutoff*SE
  ATE_CI_ub = ATE + normal.cutoff*SE
  
  # Estimate RMSE and CI coverage
  true_tau <- mean(outcome.data[, "tau"])
  rmse <- sqrt(mean((ATE - outcome.data[, "tau"])^2))
  cover <- (true_tau >= ATE_CI_lb) & (true_tau <= ATE_CI_ub)
  return(c(true_tau, ATE, ATE_CI_lb, ATE_CI_ub, rmse, cover))
}

TMLE_XBART_estimator <- function(df, type = "full", 
                                 estimate_propensity = TRUE){
  # Requires that the dataframe have the following columns:
  # Y, Z, x_1, x_2, x_3, x_4, x_5
  outcome_col <- "Y"
  treatment_col <- "Z"
  covariates <- c("x_1", "x_2", "x_3", "x_4", "x_5")
  
  if (type == "full"){
    propensity.data <- df
    outcome.data <- df
    N = length(df$Y)
  } else if (type == "ssl"){
    propensity.data <- df
    outcome.data <- df[df$M == 0, ]
    N = length(df[df$M == 0, "Y"])
  } else if (type == "complete_case"){
    propensity.data <- df[df$M == 0, ]
    outcome.data <- df[df$M == 0, ]
    N = length(df[df$M == 0, "Y"])
  }
  
  W_a = outcome.data[, c(treatment_col, covariates)]
  W_1 = cbind(Z=1, outcome.data[, covariates])
  W_0 = cbind(Z=0, outcome.data[, covariates])
  Y = outcome.data[, outcome_col]
  
  # Fit prognostic model
  get_XBART_params <- function(y, p) {
    XBART_params = list(num_trees = 30, # number of trees
                        num_sweeps = 40, # number of sweeps (samples of the forest)
                        n_min = 1, # minimal node size
                        alpha = 0.95, # BART prior parameter
                        beta = 1.25, # BART prior parameter
                        mtry = p, # number of variables sampled in each split
                        burnin = 15,
                        no_split_penality = "Auto"
    ) # burnin of MCMC sample
    num_tress = XBART_params$num_trees
    XBART_params$max_depth = 250
    XBART_params$num_cutpoints = 50;
    # number of adaptive cutpoints
    XBART_params$tau = var(y) / num_tress # prior variance of mu (leaf parameter)
    return(XBART_params)
  }
  params = get_XBART_params(Y, ncol(W_a))
  dcat = 0
  parl = F
  prognostic.model = XBART(as.matrix(Y), as.matrix(W_a), as.matrix(W_a), p_categorical = dcat,
                           params$num_trees, params$num_sweeps,
                           params$max_depth, params$n_min,
                           alpha = params$alpha, beta = params$beta,
                           tau = params$tau, s = 1, kap = 1,
                           mtry = params$mtry, verbose = verbose,
                           num_cutpoints = params$num_cutpoints,
                           parallel = parl, random_seed = 100,
                           no_split_penality = params$no_split_penality)
  Q_1_pred = rowMeans(
    predict(prognostic.model, as.matrix(W_1))[, params$burnin:params$num_sweeps]
  )
  Q_0_pred = rowMeans(
    predict(prognostic.model, as.matrix(W_0))[, params$burnin:params$num_sweeps]
  )
  
  # prognostic.model = bart(W_a, Y, keeptrees = T, verbose = F)
  # Q_a_pred = colMeans(predict(prognostic.model, W_a))
  # Q_1_pred = colMeans(predict(prognostic.model, W_1))
  # Q_0_pred = colMeans(predict(prognostic.model, W_0))
  Q_pred = cbind(Q_1_pred, Q_0_pred)
  
  # Build propensity model and predict propensity scores
  X.train = propensity.data[, covariates]
  Z.train = propensity.data[, treatment_col]
  X.pred = outcome.data[, covariates]
  Z.pred = outcome.data[, treatment_col]
  
  # Estimate tau and 95% confidence interval
  if (estimate_propensity){
    
    get_XBART_params <- function(y, p){
      XBART_params = list(M = 15,
                          L = 1,
                          nsweeps = 150,
                          Nmin = 100,
                          alpha = 0.95,
                          beta = 1.25,
                          mtry = p,
                          burnin = 15)
      num_trees = XBART_params$M
      XBART_params$max_depth = 250
      XBART_params$Ncutpoints = 100;
      XBART_params$tau = var(y)/(num_trees)
      XBART_params$a = 0.000001; XBART_params$b = 0.000001;
      return(XBART_params)
    }
    params = get_XBART_params(Z.train, ncol(X.train))
    dcat = 0
    parl = F
    
    propensity.model = XBART.Probit(
      as.matrix(Z.train), as.matrix(X.train), as.matrix(X.train), num_trees = params$M,
      L = 1, num_sweeps = params$nsweeps, 
      max_depth = params$max_depth, Nmin = 10,
      num_cutpoints = params$Ncutpoints,
      alpha = params$alpha, beta = params$beta,
      tau = params$tau, s= 1,kap = 1,
      mtry = params$mtry, p_categorical = dcat,
      draw_sigma = FALSE, m_update_sigma = TRUE,
      parallel = parl, random_seed = 10
    )
    propensity.pred = apply(
      pnorm(predict(propensity.model, as.matrix(X.pred))[,params$burnin:params$nsweeps]), 1, mean
    )
    # propensity.model = bart(X.train, Z.train, keeptrees = T, verbose = F)
    # propensity.pred = colMeans(pnorm(predict(propensity.model, X.pred)))
  } else{
    propensity.pred = outcome.data$pi_x
  }
  
  tau.model <- tmle(Y, Z.pred, X.pred, 
                    Q = Q_pred, g1W = propensity.pred)
  ATE = tau.model$estimates$ATE$psi
  ATE_CI_lb = tau.model$estimates$ATE$CI[1]
  ATE_CI_ub = tau.model$estimates$ATE$CI[2]
  
  # Estimate RMSE and CI coverage
  true_tau <- mean(outcome.data[, "tau"])
  rmse <- sqrt(mean((ATE - outcome.data[, "tau"])^2))
  cover <- (true_tau >= ATE_CI_lb) & (true_tau <= ATE_CI_ub)
  return(c(true_tau, ATE, ATE_CI_lb, ATE_CI_ub, rmse, cover))
}

BCF_XBART_estimator <- function(df, type = "full", alpha = 0.05, 
                          estimate_propensity = TRUE){
  # Requires that the dataframe have the following columns:
  # Y, Z, x_1, x_2, x_3, x_4, x_5
  outcome_col <- "Y"
  treatment_col <- "Z"
  covariates <- c("x_1", "x_2", "x_3", "x_4", "x_5")
  
  if (type == "full"){
    propensity.data <- df
    outcome.data <- df
    N = length(df$Y)
  } else if (type == "ssl"){
    propensity.data <- df
    outcome.data <- df[df$M == 0, ]
    N = length(df[df$M == 0, "Y"])
  } else if (type == "complete_case"){
    propensity.data <- df[df$M == 0, ]
    outcome.data <- df[df$M == 0, ]
    N = length(df[df$M == 0, "Y"])
  }
  
  # Build propensity model and predict propensity scores
  X.train = propensity.data[, covariates]
  Z.train = propensity.data[, treatment_col]
  X.pred = outcome.data[, covariates]
  Z.pred = outcome.data[, treatment_col]
  
  # Estimate tau and 95% confidence interval
  if (estimate_propensity){
    get_XBART_params <- function(y, p){
      XBART_params = list(M = 15,
                          L = 1,
                          nsweeps = 150,
                          Nmin = 100,
                          alpha = 0.95,
                          beta = 1.25,
                          mtry = p,
                          burnin = 15)
      num_trees = XBART_params$M
      XBART_params$max_depth = 250
      XBART_params$Ncutpoints = 100;
      XBART_params$tau = var(y)/(num_trees)
      XBART_params$a = 0.000001; XBART_params$b = 0.000001;
      return(XBART_params)
    }
    params = get_XBART_params(Z.train, ncol(X.train))
    dcat = 0
    parl = F
    
    propensity.model = XBART.Probit(
      as.matrix(Z.train), as.matrix(X.train), as.matrix(X.train), num_trees = params$M,
      L = 1, num_sweeps = params$nsweeps, 
      max_depth = params$max_depth, Nmin = 10,
      num_cutpoints = params$Ncutpoints,
      alpha = params$alpha, beta = params$beta,
      tau = params$tau, s= 1,kap = 1,
      mtry = params$mtry, p_categorical = dcat,
      draw_sigma = FALSE, m_update_sigma = TRUE,
      parallel = parl, random_seed = 10
    )
    propensity.pred = apply(
      pnorm(predict(propensity.model, as.matrix(X.pred))[,params$burnin:params$nsweeps]), 1, mean
    )
  } else{
    propensity.pred = outcome.data$pi_x
  }
  
  bcf.model <- bcf(outcome.data$Y, outcome.data$Z, as.matrix(X.pred), 
                   pihat = propensity.pred, nburn=100, nsim=100)
  ATE = mean(colMeans(bcf.model$tau))
  ATE_CI_lb = unname(quantile(colMeans(bcf.model$tau), (alpha/2)))
  ATE_CI_ub = unname(quantile(colMeans(bcf.model$tau), (1-(alpha/2))))
  
  # Estimate RMSE and CI coverage
  true_tau <- mean(outcome.data[, "tau"])
  rmse <- sqrt(mean((ATE - outcome.data[, "tau"])^2))
  cover <- (true_tau >= ATE_CI_lb) & (true_tau <= ATE_CI_ub)
  return(c(true_tau, ATE, ATE_CI_lb, ATE_CI_ub, rmse, cover))
}

XBCF_estimator <- function(df, type = "full", alpha = 0.05, 
                          estimate_propensity = TRUE){
  # Requires that the dataframe have the following columns:
  # Y, Z, x_1, x_2, x_3, x_4, x_5
  outcome_col <- "Y"
  treatment_col <- "Z"
  covariates <- c("x_1", "x_2", "x_3", "x_4", "x_5")
  
  if (type == "full"){
    propensity.data <- df
    outcome.data <- df
    N = length(df$Y)
  } else if (type == "ssl"){
    propensity.data <- df
    outcome.data <- df[df$M == 0, ]
    N = length(df[df$M == 0, "Y"])
  } else if (type == "complete_case"){
    propensity.data <- df[df$M == 0, ]
    outcome.data <- df[df$M == 0, ]
    N = length(df[df$M == 0, "Y"])
  }
  
  # Build propensity model and predict propensity scores
  X.train = propensity.data[, covariates]
  Z.train = propensity.data[, treatment_col]
  X.pred = outcome.data[, covariates]
  Z.pred = outcome.data[, treatment_col]
  
  # Estimate tau and 95% confidence interval
  if (estimate_propensity){
    get_XBART_params <- function(y, p){
      XBART_params = list(M = 15,
                          L = 1,
                          nsweeps = 150,
                          Nmin = 100,
                          alpha = 0.95,
                          beta = 1.25,
                          mtry = p,
                          burnin = 15)
      num_trees = XBART_params$M
      XBART_params$max_depth = 250
      XBART_params$Ncutpoints = 100;
      XBART_params$tau = var(y)/(num_trees)
      XBART_params$a = 0.000001; XBART_params$b = 0.000001;
      return(XBART_params)
    }
    params = get_XBART_params(Z.train, ncol(X.train))
    dcat = 0
    parl = F
    
    propensity.model = XBART.Probit(
      as.matrix(Z.train), as.matrix(X.train), as.matrix(X.train), num_trees = params$M,
      L = 1, num_sweeps = params$nsweeps, 
      max_depth = params$max_depth, Nmin = 10,
      num_cutpoints = params$Ncutpoints,
      alpha = params$alpha, beta = params$beta,
      tau = params$tau, s= 1,kap = 1,
      mtry = params$mtry, p_categorical = dcat,
      draw_sigma = FALSE, m_update_sigma = TRUE,
      parallel = parl, random_seed = 10
    )
    propensity.pred = apply(
      pnorm(predict(propensity.model, as.matrix(X.pred))[,params$burnin:params$nsweeps]), 1, mean
    )
    # propensity.model = bart(X.train, Z.train, keeptrees = T, verbose = F)
    # propensity.pred = colMeans(pnorm(predict(propensity.model, X.pred)))
  } else{
    propensity.pred = outcome.data$pi_x
  }
  
  # Prepare XBCF run
  Y = df$Y - mean(df$Y)
  sdy = sd(Y)
  Y = Y/sdy
  
  burnin = 25
  sweeps = 100
  treesmu = 60
  treestau = 30
  
  tau1 = 0.9*var(Y)/treesmu
  tau2 = 0.1*var(Y)/treestau
  
  xbcf_fit = XBCF(
    as.matrix(Y), as.matrix(X.pred), as.matrix(X.pred), as.matrix(Z.pred), num_sweeps = sweeps, 
    burnin = burnin, max_depth = 50, Nmin = 1, num_cutpoints = 30, no_split_penality = "Auto",
    mtry_pr = ncol(X.pred), mtry_trt = ncol(X.pred), p_categorical_pr = 5,  p_categorical_trt = 5,
    num_trees_pr = treesmu, alpha_pr = 0.95, beta_pr = 1.25, tau_pr = tau1, kap_pr = 1, s_pr = 1, 
    pr_scale = FALSE, num_trees_trt = treestau, alpha_trt = 0.25, beta_trt = 2, tau_trt = tau2, 
    kap_trt =1, s_trt = 1, trt_scale = FALSE, verbose = FALSE, a_scaling = TRUE, b_scaling = TRUE
  )
  # compute tauhats as (b1-b0)*tau
  th = xbcf_fit$tauhats
  b = xbcf_fit$b_draws
  seq <- (burnin+1):sweeps
  for (i in seq){ th[,i] = th[,i] * (b[i,2] - b[i,1]) }
  tauhats = rowSums(th[,(burnin+1):sweeps])/(sweeps-burnin)
  tauhats = tauhats*sdy
  ATE = mean(tauhats)
  ATE_CI_lb = unname(quantile(tauhats, (alpha/2)))
  ATE_CI_ub = unname(quantile(tauhats, (1-(alpha/2))))
  
  # bcf.model <- bcf(outcome.data$Y, outcome.data$Z, as.matrix(X.pred), 
  #                  pihat = propensity.pred, nburn=100, nsim=100)
  # ATE = mean(colMeans(bcf.model$tau))
  # ATE_CI_lb = unname(quantile(colMeans(bcf.model$tau), (alpha/2)))
  # ATE_CI_ub = unname(quantile(colMeans(bcf.model$tau), (1-(alpha/2))))
  # 
  # Estimate RMSE and CI coverage
  true_tau <- mean(outcome.data[, "tau"])
  rmse <- sqrt(mean((tauhats - outcome.data[, "tau"])^2))
  cover <- (true_tau >= ATE_CI_lb) & (true_tau <= ATE_CI_ub)
  return(c(true_tau, ATE, ATE_CI_lb, ATE_CI_ub, rmse, cover))
}