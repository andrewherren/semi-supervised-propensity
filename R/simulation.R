## Functions to generate simulated data

bcf_simulations <- function(n = 250, linear = TRUE, homogeneous = TRUE,
                            missing_mechanism = c("None", "MCAR", 
                                                  "linear", "truncated"), 
                            missing_pct = 0.5, confounded = TRUE){
  # covariates
  x_1 = rnorm(n, 0, 1)
  x_2 = rnorm(n, 0, 1)
  x_3 = rnorm(n, 0, 1)
  x_4 = rbinom(n, 1, 0.5)
  x_5 = apply((rmultinom(n, 1, c(0.25, 0.5, 0.25))==1), 2, which)
  
  # transformation of x_5
  g <- function(x_5){
    ifelse(x_5 == 1, 2, ifelse(x_5 == 2, -1, -4))
  }
  
  # prognostic function(s)
  if (linear){
    mu = 1 + g(x_5) + x_1*x_3
  } else{
    mu = -6 + g(x_5) + 6*abs(x_3-1)
  }
  
  # propensity function
  if (confounded){
    s = sd(mu)
    u = runif(n, 0, 1)
    pi_x = 0.8*(pnorm((3*mu/s) - 0.5*x_1) + 0.05 + u/10)
  } else{
    pi_x = 0.5
  }

  # treatment effect
  if (homogeneous){
    tau = 3
  } else{
    tau = 1 + 2*x_2*x_4
  }
  
  # generate outcome
  z = rbinom(n, 1, pi_x)
  y = mu + tau*z + rnorm(n, 0, 1)
  
  # generate missing values
  if (missing_mechanism=="None") {
    M = rep(0, n)
  } else if (missing_mechanism=="MCAR"){
    # M = rbinom(n, 1, missing_pct)
    n_missing = floor((missing_pct)*n)
    M = sample(c(rep(1,n_missing),rep(0,n-n_missing)))
  } else if (missing_mechanism=="linear"){
    # Designed to make missingness correlated with x_2
    # and have a total % unlabeled of roughly 90%
    q = 1/(1 + exp(-(2.40 + 0.7*x_2)))
    M = rbinom(n, 1, q)
  } else if (missing_mechanism=="truncated"){
    # # Taking the center 10% of a standard normal distribution
    # M = ifelse(abs(x_2) <= 0.1256613, 0, 1)
    # Exclude the outermost missing_pct% of Y values
    lower_bound = quantile(y, (missing_pct/2))
    upper_bound = quantile(y, 1 - (missing_pct/2))
    M = ifelse((y >= lower_bound) & (y <= upper_bound), 0, 1)
  } else {
    M = rep(0, n)
  }
  
  return(list(Y = y, X = cbind(x_1, x_2, x_3, x_4, x_5), 
              Z = z, M = M, tau = tau, pi_x = pi_x))
}
