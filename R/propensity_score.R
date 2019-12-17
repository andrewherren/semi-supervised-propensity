## Estimating propensity scores using a variety of predictive models

propensity.estimates <- function(treatment, X, X.treatment=X){
  # Estimate propensity scores using a variety of methods
  # @param treatment: vector of treatment assignments
  # @param X: vector of simulated covariates
  # @param X.treatment: vector of covariates used to assign treatment
  # @return: list of estimated treatment probabilities
  
  # Logistic regression
  # Correctly specified using only the variables predictive of tau
  correct.logit.W <- glm(treatment ~ X.treatment, 
                         family = binomial(link = "logit"))
  W.hat.logit.correct <- predict(correct.logit.W, type="response")
  # Incorrectly specified, using all variables
  kitchen.sink.logit.W <- glm(treatment ~ X, 
                              family = binomial(link = "logit"))
  W.hat.logit.kitchen.sink <- predict(kitchen.sink.logit.W, type="response")
  
  # LASSO regression
  lasso.W <- glmnet(x = X, y = treatment, family = "binomial", alpha = 1)
  W.hat.lasso <- apply(predict(lasso.W, newx=X, type="response"), 1, mean)
  
  # Ridge regression
  ridge.W <- glmnet(x = X, y = treatment, family = "binomial", alpha = 0)
  W.hat.ridge <- apply(predict(ridge.W, newx=X, type="response"), 1, mean)
  
  # GRF - generalized random forest
  forest.W <- regression_forest(X, treatment, tune.parameters = TRUE)
  W.hat.grf <- predict(forest.W)$predictions
  
  # BART - Bayesian additive regression trees
  # bart.W <- pbart(x.train = X, y.train = treatment, 
  #                 ndpost = 5000, nskip = 2000)
  # W.hat.bart <- predict(bart.W, newdata=X)$prob.test.mean
  
  # XBART
  # print(length(treatment))
  # print(length(X))
  # XBART.params = get_XBART_params(treatment)
  # XBART.w = XBART.Probit(y = as.matrix(treatment), X = as.matrix(X), 
  #                        Xtest = as.matrix(X), p_categorical = 0, 
  #                        XBART.params$num_trees, XBART.params$num_sweeps, 
  #                        XBART.params$max_depth, XBART.params$n_min, 
  #                        alpha = XBART.params$alpha, beta = XBART.params$beta, 
  #                        tau = XBART.params$tau, s = 1, kap = 1, 
  #                        mtry = ncol(as.matrix(X))-1, verbose = verbose, 
  #                        num_cutpoints = XBART.params$num_cutpoints, 
  #                        parallel = parl, random_seed = 100, 
  #                        no_split_penality = XBART.params$no_split_penality)
  # W.hat.xbart <- XBART.w$yhats_test
  
  return(list(
    W.hat.logit.correct=W.hat.logit.correct, 
    W.hat.logit.kitchen.sink=W.hat.logit.kitchen.sink, 
    W.hat.lasso=W.hat.lasso, 
    W.hat.ridge=W.hat.ridge, 
    W.hat.grf=W.hat.grf 
    # W.hat.bart=W.hat.bart, 
    # W.hat.xbart=W.hat.xbart
  ))
}

mse <- function(y.hat, y){
  # Estimate mean squared error
  # @param y.hat: fitted values from a model
  # @param y: true values
  # @return: point estimate of MSE
  return(sum((y-y.hat)^2)/length(y.hat))
}

evaluate_propensities <- function(propensity.list, true.prob){
  # Evaluate propensity score estimates provided in list against
  # simulate true probability of receiving treatment
  # @param propensity_list: list of estimate treatment probabilities
  # @param true.prob: simulated true probability of receiving treatment
  # @return: mean squared error for each method provided in propensity.list
  return(sapply(propensity.list, mse, y=true.prob))
}
