########################################
# BCF simulation
########################################

# load packages
pkg.list <- c("bcf", "dbarts", "ranger", "foreach", "doParallel")
lapply(pkg.list, require, character.only = TRUE)

# set random seed
set.seed(1234)

# Set up working environment
require(here)
project_dir <- here()

# Register parallel backend
numCores <- detectCores()
registerDoParallel(numCores)

# Create new "snapshots" folder, if doesn't exist
datestamp = format(Sys.time(), "%Y%m%d%H%M")
logfile_name = file.path(project_dir, "log", paste0("R_batch_log_", datestamp, ".log"))
output_file = file(logfile_name, "w")
snapshots_subfolder = file.path(project_dir, "outputs", "snapshots", datestamp)
ifelse(!dir.exists(snapshots_subfolder), dir.create(snapshots_subfolder), FALSE)

# data size, treatment effect, and epison_y
n = 500
p_y = 2
p_z = 1
p_yz = 4
tau = 3
sigma_y = 1

# simulation parameters
beta_x_z_levels = c(0, 1, 10); n_beta_x_z <- length(beta_x_z_levels)
n_sim = 2
sim_results <- as.data.frame(matrix(rep(0, n_beta_x_z*n_sim*24), ncol = 24))
colnames(sim_results) <- c("tau_hat_case_1", "tau_hat_case_2", "tau_hat_case_3", "tau_hat_case_4", 
                           "std_tau_hat_case_1", "std_tau_hat_case_2", "std_tau_hat_case_3", "std_tau_hat_case_4", 
                           "tau_coverage_case_1", "tau_coverage_case_2", "tau_coverage_case_3", "tau_coverage_case_4", 
                           "tau_bias_case_1", "tau_bias_case_2", "tau_bias_case_3", "tau_bias_case_4", 
                           "tau_rmse_case_1", "tau_rmse_case_2", "tau_rmse_case_3", "tau_rmse_case_4", 
                           "cor_x_z_y_incl", "cor_x_z_y_excl", "cor_x_z_z", "beta")

for (b_x_z in 1:n_beta_x_z){
  t = proc.time()[3]
  
  beta_sim_result <- foreach (sim=1:n_sim, .combine = rbind) %dopar% {
    ### Covariates
    X_y <- matrix(rnorm(n*p_y, 0, 1), ncol = p_y)
    X_z <- matrix(rnorm(n*p_z, 0, 1), ncol = p_z)
    X_yz <- matrix(rnorm(n*p_yz, 0, 1), ncol = p_yz)
    colnames(X_y) <- paste0("x", 1:p_y)
    colnames(X_z) <- paste0("x", (p_y+1):(p_y+p_z))
    colnames(X_yz) <- paste0("x", (p_y+p_z+1):(p_y+p_z+p_yz))
    
    ### Treatment
    gamma <- runif(p_yz, -0.3, 0.3)
    gamma_x_z <- beta_x_z_levels[b_x_z]
    true_propensity <- pnorm(X_yz%*%gamma + X_z%*%gamma_x_z)
    # z <- (true_propensity>0.5)*1
    z <- as.matrix(rbinom(n, 1, true_propensity))
    colnames(z) <- "z"
    
    ### Outcome
    beta_excl_x_z <- runif(p_y + p_yz, -0.5, 0.5)
    beta_x_z <- 0.5
    eps <- rnorm(n, 0, sigma_y)
    y_0 <- cbind(X_y, X_yz)%*%beta_excl_x_z
    y_1 <- y_0 + tau*z
    y_excl_x_z <- z*y_1 + (1-z)*y_0 + eps
    y_incl_x_z <- y_excl_x_z + X_z%*%beta_x_z
    colnames(y_excl_x_z) <- "y"
    colnames(y_incl_x_z) <- "y"

    ### Propensity model
    # Case 1: x_z included in pihat
    # propensity.model.all <- ranger(z ~ ., data = as.data.frame(cbind(z, X_y, X_z, X_yz)))
    # propensities.all <- predict(propensity.model.all, data = as.data.frame(cbind(X_y, X_z, X_yz)), type = "response")$predictions
    propensity.model.all <- bart(cbind(X_y, X_z, X_yz), z, keeptrees = T)
    propensities.all <- colMeans(pnorm(predict(propensity.model.all, cbind(X_y, X_z, X_yz))))

    # Case 2: x_z not included in pihat
    # propensity.model.excl.x.z <- ranger(z ~ ., data = as.data.frame(cbind(z, X_y, X_yz)))
    # propensities.excl.x.z <- predict(propensity.model.excl.x.z, data = as.data.frame(cbind(X_y, X_yz)), type = "response")$predictions
    propensity.model.excl.x.z <- bart(cbind(X_y, X_yz), z, keeptrees = T)
    propensities.excl.x.z <- colMeans(pnorm(predict(propensity.model.excl.x.z, cbind(X_y, X_yz))))

    ### Treatment effect estimation
    # Case 1: x_z is confounder, included in pihat
    bcf_obj_1 = bcf(y=y_incl_x_z, z=z, x_control = cbind(X_y, X_z, X_yz), pihat=propensities.all, nburn = 100, nsim = 100, nthin = 10)
    tau.hat.bcf.1 = mean(colMeans(bcf_obj_1$tau))
    std.tau.hat.bcf.1 = sd(colMeans(bcf_obj_1$tau))
    
    # Case 2: x_j not confounder, included in pihat
    bcf_obj_2 = bcf(y=y_excl_x_z, z=z, x_control = cbind(X_y, X_z, X_yz), pihat=propensities.all, nburn = 100, nsim = 100, nthin = 10)
    tau.hat.bcf.2 = mean(colMeans(bcf_obj_2$tau))
    std.tau.hat.bcf.2 = sd(colMeans(bcf_obj_2$tau))
    
    # Case 3: x_j is confounder, not included in pihat
    bcf_obj_3 = bcf(y=y_incl_x_z, z=z, x_control = cbind(X_y, X_z, X_yz), pihat=propensities.excl.x.z, nburn = 100, nsim = 100, nthin = 10)
    tau.hat.bcf.3 = mean(colMeans(bcf_obj_3$tau))
    std.tau.hat.bcf.3 = sd(colMeans(bcf_obj_3$tau))
    
    # Case 4: x_j not confounder, not included in pihat
    bcf_obj_4 = bcf(y=y_excl_x_z, z=z, x_control = cbind(X_y, X_z, X_yz), pihat=propensities.excl.x.z, nburn = 100, nsim = 100, nthin = 10)
    tau.hat.bcf.4 = mean(colMeans(bcf_obj_4$tau))
    std.tau.hat.bcf.4 = sd(colMeans(bcf_obj_4$tau))

    # Interval coverage
    tau.coverage.1 <- ((tau >= quantile(rowMeans(bcf_obj_1$tau), 0.05)) & (tau <= quantile(rowMeans(bcf_obj_1$tau), 0.95)))*1
    tau.coverage.2 <- ((tau >= quantile(rowMeans(bcf_obj_2$tau), 0.05)) & (tau <= quantile(rowMeans(bcf_obj_2$tau), 0.95)))*1
    tau.coverage.3 <- ((tau >= quantile(rowMeans(bcf_obj_3$tau), 0.05)) & (tau <= quantile(rowMeans(bcf_obj_3$tau), 0.95)))*1
    tau.coverage.4 <- ((tau >= quantile(rowMeans(bcf_obj_4$tau), 0.05)) & (tau <= quantile(rowMeans(bcf_obj_4$tau), 0.95)))*1
    
    # Bias
    tau.bias.1 <- tau.hat.bcf.1-tau
    tau.bias.2 <- tau.hat.bcf.2-tau
    tau.bias.3 <- tau.hat.bcf.3-tau
    tau.bias.4 <- tau.hat.bcf.4-tau
    
    # Bias
    tau.rmse.1 <- sqrt((tau.hat.bcf.1-tau)^2)
    tau.rmse.2 <- sqrt((tau.hat.bcf.2-tau)^2)
    tau.rmse.3 <- sqrt((tau.hat.bcf.3-tau)^2)
    tau.rmse.4 <- sqrt((tau.hat.bcf.4-tau)^2)
    
    # Compile all simulation results
    data.frame(tau_hat_case_1 = tau.hat.bcf.1, tau_hat_case_2 = tau.hat.bcf.2, tau_hat_case_3 = tau.hat.bcf.3, tau_hat_case_4 = tau.hat.bcf.4, 
         std_tau_hat_case_1 = std.tau.hat.bcf.1, std_tau_hat_case_2 = std.tau.hat.bcf.2, std_tau_hat_case_3 = std.tau.hat.bcf.3, std_tau_hat_case_4 = std.tau.hat.bcf.4, 
         tau_coverage_case_1 = tau.coverage.1, tau_coverage_case_2 = tau.coverage.2, tau_coverage_case_3 = tau.coverage.3, tau_coverage_case_4 = tau.coverage.4, 
         tau_bias_case_1 = tau.bias.1, tau_bias_case_2 = tau.bias.2, tau_bias_case_3 = tau.bias.3, tau_bias_case_4 = tau.bias.4, 
         tau_rmse_case_1 = tau.rmse.1, tau_rmse_case_2 = tau.rmse.2, tau_rmse_case_3 = tau.rmse.3, tau_rmse_case_4 = tau.rmse.4, 
         cor_x_z_y_incl = cor(y_incl_x_z, X_z), cor_x_z_y_excl = cor(y_excl_x_z, X_z), cor_x_z_z = cor(z, X_z), beta = beta_x_z_levels[b_x_z])
  }
  sim_results[((b_x_z-1)*n_sim + 1):(b_x_z*n_sim), ] <- beta_sim_result
  
  t = proc.time()[3] - t
  cat("==================================================\n", file = output_file)
  cat(paste("Run time:", t, "\n"), file = output_file)
  cat(paste("Iteration", b_x_z, "of", n_beta_x_z, "complete\n\n"), file = output_file)
}
close(output_file)

# Save outputs to file
write.csv(sim_results, file.path(project_dir, "outputs", "snapshots", datestamp, paste0("BCF_sim_", datestamp, ".csv")))
