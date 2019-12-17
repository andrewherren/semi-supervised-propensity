project_dir <- "~/Github/RTG-2019-2H"
setwd(project_dir)
source("R/simulation.R", local=FALSE, echo=TRUE)
set.seed(1234)
n_sim = 250

save_simulations <- function(n = 5000, linear = TRUE, 
                             homogeneous = TRUE, 
                             missing_mechanism = "MCAR", 
                             missing_pct = 0.90, 
                             confounded = TRUE){
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
  
  # Save to CSV
  lin = ifelse(linear, "lin", "nonlin")
  hom = ifelse(homogeneous, "hom", "het"); sample = n
  missing = paste(missing_mechanism, floor(missing_pct*100), sep = "_")
  conf = ifelse(confounded, "conf", "unconf")
  subfolder_name = paste(lin, hom, sample, missing, conf, sep="_")
  file_name = paste0("df", i, ".csv")
  
  # Save down
  directory_path = file.path(project_dir, "data", "simulations", subfolder_name)
  # Create sub-directory if doesn't already exist
  ifelse(!dir.exists(directory_path), dir.create(directory_path), FALSE)
  # Save simulated dataset to CSV
  data_file <- file.path(directory_path, file_name)
  write.csv(x = causal_data, file = data_file, row.names = F)
}

for (i in 1:n_sim){
  ###########################################
  # CASE 1
  n = 5000; linear = TRUE; homogeneous = TRUE
  missing_mechanism = "truncated"; missing_pct = 0.90
  confounded = TRUE
  save_simulations(n = n, linear = linear, 
                   homogeneous = homogeneous, 
                   missing_mechanism = missing_mechanism, 
                   missing_pct = missing_pct, 
                   confounded = confounded)
  
  ###########################################
  # CASE 2
  n = 5000; linear = TRUE; homogeneous = TRUE
  missing_mechanism = "truncated"; missing_pct = 0.90
  confounded = FALSE
  save_simulations(n = n, linear = linear, 
                   homogeneous = homogeneous, 
                   missing_mechanism = missing_mechanism, 
                   missing_pct = missing_pct, 
                   confounded = confounded)
  
  ###########################################
  # CASE 3
  n = 1000; linear = TRUE; homogeneous = TRUE
  missing_mechanism = "truncated"; missing_pct = 0.90
  confounded = TRUE
  save_simulations(n = n, linear = linear, 
                   homogeneous = homogeneous, 
                   missing_mechanism = missing_mechanism, 
                   missing_pct = missing_pct, 
                   confounded = confounded)
  
  ###########################################
  # CASE 4
  n = 1000; linear = TRUE; homogeneous = TRUE
  missing_mechanism = "truncated"; missing_pct = 0.90
  confounded = FALSE
  save_simulations(n = n, linear = linear, 
                   homogeneous = homogeneous, 
                   missing_mechanism = missing_mechanism, 
                   missing_pct = missing_pct, 
                   confounded = confounded)
  
  ###########################################
  # CASE 5
  n = 500; linear = TRUE; homogeneous = TRUE
  missing_mechanism = "truncated"; missing_pct = 0.90
  confounded = TRUE
  save_simulations(n = n, linear = linear, 
                   homogeneous = homogeneous, 
                   missing_mechanism = missing_mechanism, 
                   missing_pct = missing_pct, 
                   confounded = confounded)
  
  ###########################################
  # CASE 6
  n = 500; linear = TRUE; homogeneous = TRUE
  missing_mechanism = "truncated"; missing_pct = 0.90
  confounded = FALSE
  save_simulations(n = n, linear = linear, 
                   homogeneous = homogeneous, 
                   missing_mechanism = missing_mechanism, 
                   missing_pct = missing_pct, 
                   confounded = confounded)
  
  ###########################################
  # CASE 7
  n = 250; linear = TRUE; homogeneous = TRUE
  missing_mechanism = "truncated"; missing_pct = 0.90
  confounded = TRUE
  save_simulations(n = n, linear = linear, 
                   homogeneous = homogeneous, 
                   missing_mechanism = missing_mechanism, 
                   missing_pct = missing_pct, 
                   confounded = confounded)
  
  ###########################################
  # CASE 8
  n = 250; linear = TRUE; homogeneous = TRUE
  missing_mechanism = "truncated"; missing_pct = 0.90
  confounded = FALSE
  save_simulations(n = n, linear = linear, 
                   homogeneous = homogeneous, 
                   missing_mechanism = missing_mechanism, 
                   missing_pct = missing_pct, 
                   confounded = confounded)
}