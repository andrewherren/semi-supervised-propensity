# Compare approaches on NHEFS data

# Load necessary packages
pkg.list <- c("dbarts", "bcf", "here")
lapply(pkg.list, require, character.only = TRUE)

# Set random seed
set.seed(2020)

# Set up access to project directory
project_dir <- here()
datestamp = format(Sys.time(), "%Y%m%d%H%M")
logfile_name = file.path(project_dir, "log", paste0("R_batch_log_", datestamp, ".log"))
output_file = file(logfile_name, "w")
nhefs_data_file = file.path(project_dir, "data", "raw", "nhefs.csv")
nhefs_data_url = "https://cdn1.sph.harvard.edu/wp-content/uploads/sites/1268/1268/20/nhefs.csv"

load_nhefs_data <- function(csv_file_name, data_url){
  if (!file.exists(csv_file_name)) {
    # Download file if not already on disk
    download.file(data_url, destfile = csv_file_name)
  }
  
  # Read from downloaded CSV file
  nhefs_data <- read.csv(csv_file_name, header = TRUE)
  return(nhefs_data)
}

IPW_BART_NHEFS <- function(nhefs_data, M = -1, alpha = 0.05, type = "full"){
  
  # Confounders used in chapter 12 of Hernan and Robins (2020)
  confounders <- c("sex", "race", "age", "smokeintensity", "smokeyrs",
                   "wt71", "education", "exercise", "active")
  # Treatment variable: smoking cessation
  treatment <- "qsmk"
  # Outcome of interest: weight gain between 1971 and 1982
  outcome <- "wt82_71"
  analysis_cols <- c(confounders, treatment, outcome)
  
  if (type == "full"){
    propensity.data <- nhefs_data[, analysis_cols]
    outcome.data <- nhefs_data[, analysis_cols]
    N = nrow(nhefs_data)
  } else if (type == "ssl"){
    stopifnot(M != -1)
    propensity.data <- nhefs_data[, analysis_cols]
    outcome.data <- nhefs_data[M == 0, analysis_cols]
    N = nrow(nhefs_data[M == 0, ])
  } else if (type == "complete_case"){
    stopifnot(M != -1)
    propensity.data <- nhefs_data[M == 0, analysis_cols]
    outcome.data <- nhefs_data[M == 0, analysis_cols]
    N = nrow(nhefs_data[M == 0, ])
  }
  
  # Fit propensity model
  propensity_model <- bart(propensity.data[, confounders], 
                           propensity.data[, treatment], 
                           keeptrees = T, verbose = F)
  # Estimate propensity scores on the outcome data
  propensity_scores = colMeans(pnorm(predict(propensity_model, 
                                             outcome.data[, confounders])))

  # Estimate ATE
  X = outcome.data[, confounders]
  Z = outcome.data[, treatment]
  Y = outcome.data[, outcome]
  ATE.summand = ((Z - propensity_scores)*Y)/
    (propensity_scores*(1 - propensity_scores))
  ATE = mean(ATE.summand)
  
  # Estimate variance of ATE
  score = (Z - propensity_scores)*X
  ATE.df = data.frame(ATE.summand = ATE.summand, score = score)
  ATE.score.model <- lm(ATE.summand ~ ., data=ATE.df)
  SE = sqrt(sum((ATE.score.model$residuals)^2)/N)/sqrt(N)
  
  # Return estimand and confidence interval
  normal.cutoff <- qnorm(1-(alpha/2))
  ATE_CI_lb = ATE - normal.cutoff*SE
  ATE_CI_ub = ATE + normal.cutoff*SE
  return(c(ATE_CI_lb, ATE, ATE_CI_ub))
}

BCF_NHEFS <- function(nhefs_data, M = -1, alpha = 0.05, type = "full"){
  
  # Confounders used in chapter 12 of Hernan and Robins (2020)
  confounders <- c("sex", "race", "age", "smokeintensity", "smokeyrs",
                   "wt71", "education", "exercise", "active")
  # Treatment variable: smoking cessation
  treatment <- "qsmk"
  # Outcome of interest: weight gain between 1971 and 1982
  outcome <- "wt82_71"
  analysis_cols <- c(confounders, treatment, outcome)
  
  if (type == "full"){
    propensity.data <- nhefs_data[, analysis_cols]
    outcome.data <- nhefs_data[, analysis_cols]
    N = nrow(nhefs_data)
  } else if (type == "ssl"){
    stopifnot(M != -1)
    propensity.data <- nhefs_data[, analysis_cols]
    outcome.data <- nhefs_data[M == 0, analysis_cols]
    N = nrow(nhefs_data[M == 0, ])
  } else if (type == "complete_case"){
    stopifnot(M != -1)
    propensity.data <- nhefs_data[M == 0, analysis_cols]
    outcome.data <- nhefs_data[M == 0, analysis_cols]
    N = nrow(nhefs_data[M == 0, ])
  }
  
  # Fit propensity model
  propensity_model <- bart(propensity.data[, confounders], 
                           propensity.data[, treatment], 
                           keeptrees = T, verbose = F)
  # Estimate propensity scores on the outcome data
  propensity_scores = colMeans(pnorm(predict(propensity_model, 
                                             outcome.data[, confounders])))

  # Fit BCF model
  X = outcome.data[, confounders]
  Z = outcome.data[, treatment]
  Y = outcome.data[, outcome]
  bcf.model <- bcf(Y, Z, as.matrix(X), pihat = propensity_scores, 
                   nburn=100, nsim=100)
  
  # Estimate ATE and variance of ATE
  ATE = mean(colMeans(bcf.model$tau))
  ATE_CI_lb = unname(quantile(colMeans(bcf.model$tau), (alpha/2)))
  ATE_CI_ub = unname(quantile(colMeans(bcf.model$tau), (1-(alpha/2))))
  return(c(ATE_CI_lb, ATE, ATE_CI_ub))
}

nhefs_data_experiment <- function(nhefs_data, missing_pct, alpha = 0.05){
  # Select some data to be marked as missing
  nhefs.reduced.sample.size <- nrow(nhefs_data)
  missing <- rbinom(nhefs.reduced.sample.size, 1, missing_pct)
  
  # Run the complete case approach using both IPW-Bart and BCF
  IPW_BART_CC <- IPW_BART_NHEFS(nhefs_data, M = missing, alpha = alpha, 
                                type = "complete_case")
  BCF_CC <- BCF_NHEFS(nhefs_data, M = missing, alpha = alpha, 
                      type = "complete_case")
  
  # Run the semi-supervised approach using both IPW-Bart and BCF
  IPW_BART_SSL <- IPW_BART_NHEFS(nhefs_data, M = missing, alpha = alpha, 
                                 type = "ssl")
  BCF_SSL <- BCF_NHEFS(nhefs_data, M = missing, alpha = alpha, 
                       type = "ssl")
  
  # Return results in tabular format
  default_col_labels <- c("CI_lb", "ATE", "CI_ub")
  IPW_BART_CC_labels <- paste(default_col_labels, "_IPW_BART_CC", sep = "")
  IPW_BART_SSL_labels <- paste(default_col_labels, "_IPW_BART_SSL", sep = "")
  BCF_CC_labels <- paste(default_col_labels, "_BCF_CC", sep = "")
  BCF_SSL_labels <- paste(default_col_labels, "_BCF_SSL", sep = "")
  output_vector <- c(IPW_BART_CC, IPW_BART_SSL, BCF_CC, BCF_SSL)
  names(output_vector) <- c(IPW_BART_CC_labels, IPW_BART_SSL_labels, 
                            BCF_CC_labels, BCF_SSL_labels)
  return(output_vector)
}

##############################################################
# Data experiment: ATE of smoking cessation on weight gain
# using subset of NHEFS and deleting outcomes at random
# to test the semi-supervised approach

# Import NHEFS data
nhefs.raw <- load_nhefs_data(nhefs_data_file, nhefs_data_url)

# Prepare data for analysis
# Restrict to individuals with non-missing values of
# sex, age, race, height, weight, education, alcohol use,
# and smoking intensity in both years (71 and 82)
nhefs <- nhefs.raw[
  (!is.na(nhefs.raw$age)) &
  (!is.na(nhefs.raw$sex)) &
  (!is.na(nhefs.raw$race)) &
  (!is.na(nhefs.raw$wt71)) &
  (!is.na(nhefs.raw$wt82)) &
  (!is.na(nhefs.raw$ht)) &
  (!is.na(nhefs.raw$education)) &
  (!is.na(nhefs.raw$alcoholfreq)) &
  (!is.na(nhefs.raw$smokeintensity)) &
  (!is.na(nhefs.raw$smkintensity82_71))
  ,
]
# Check that averages match those reported in book
#colMeans(nhefs[nhefs$qsmk == 1, c('age', 'sex', 'race', 'wt71', 'smokeintensity')])
#colMeans(nhefs[nhefs$qsmk == 0, c('age', 'sex', 'race', 'wt71', 'smokeintensity')])

# Estimate the ATE on the full data using both BCF and IPW-Bart
# Hernan and Robins reports an estimated ATE of 3.4 and 
# 95% confidence interval of 2.4 - 4.5
IPW_BART_FULL_DATA <- IPW_BART_NHEFS(nhefs, M = missing, alpha = 0.05, 
                                     type = "full")
BCF_FULL_DATA <- BCF_NHEFS(nhefs, M = missing, alpha = 0.05, 
                           type = "full")

# Run repeated iterations of the complete case and semi-supervised
# approached with different levels of data missing 
n_replicates <- 200
missing_80_replicate_study <- replicate(n_replicates, 
  nhefs_data_experiment(nhefs, missing_pct = 0.8, alpha = 0.05))
missing_50_replicate_study <- replicate(n_replicates, 
  nhefs_data_experiment(nhefs, missing_pct = 0.5, alpha = 0.05))
missing_20_replicate_study <- replicate(n_replicates, 
  nhefs_data_experiment(nhefs, missing_pct = 0.2, alpha = 0.05))

# Compute average ATEs for each method across each missing percentage
BCF_CC_ATE = c(ATE_80 = mean(missing_80_replicate_study["ATE_BCF_CC", ]), 
               ATE_50 = mean(missing_50_replicate_study["ATE_BCF_CC", ]), 
               ATE_20 = mean(missing_20_replicate_study["ATE_BCF_CC", ]))
BCF_SSL_ATE = c(ATE_80 = mean(missing_80_replicate_study["ATE_BCF_SSL", ]), 
                ATE_50 = mean(missing_50_replicate_study["ATE_BCF_SSL", ]), 
                ATE_20 = mean(missing_20_replicate_study["ATE_BCF_SSL", ]))
IPW_BART_CC_ATE = c(ATE_80 = mean(missing_80_replicate_study["ATE_IPW_BART_CC", ]), 
                    ATE_50 = mean(missing_50_replicate_study["ATE_IPW_BART_CC", ]), 
                    ATE_20 = mean(missing_20_replicate_study["ATE_IPW_BART_CC", ]))
IPW_BART_SSL_ATE = c(ATE_80 = mean(missing_80_replicate_study["ATE_IPW_BART_SSL", ]), 
                     ATE_50 = mean(missing_50_replicate_study["ATE_IPW_BART_SSL", ]), 
                     ATE_20 = mean(missing_20_replicate_study["ATE_IPW_BART_SSL", ]))

# Compute average interval coverage for each method across each missing percentage
BCF_CC_cover = c(cover_80 = mean((missing_80_replicate_study["CI_lb_BCF_CC", ] < BCF_FULL_DATA[2]) & 
                                   (missing_80_replicate_study["CI_ub_BCF_CC", ] > BCF_FULL_DATA[2])), 
                 cover_50 = mean((missing_50_replicate_study["CI_lb_BCF_CC", ] < BCF_FULL_DATA[2]) & 
                                   (missing_50_replicate_study["CI_ub_BCF_CC", ] > BCF_FULL_DATA[2])), 
                 cover_20 = mean((missing_20_replicate_study["CI_lb_BCF_CC", ] < BCF_FULL_DATA[2]) & 
                                   (missing_20_replicate_study["CI_ub_BCF_CC", ] > BCF_FULL_DATA[2])))
BCF_SSL_cover = c(cover_80 = mean((missing_80_replicate_study["CI_lb_BCF_SSL", ] < BCF_FULL_DATA[2]) & 
                                   (missing_80_replicate_study["CI_ub_BCF_SSL", ] > BCF_FULL_DATA[2])), 
                 cover_50 = mean((missing_50_replicate_study["CI_lb_BCF_SSL", ] < BCF_FULL_DATA[2]) & 
                                   (missing_50_replicate_study["CI_ub_BCF_SSL", ] > BCF_FULL_DATA[2])), 
                 cover_20 = mean((missing_20_replicate_study["CI_lb_BCF_SSL", ] < BCF_FULL_DATA[2]) & 
                                   (missing_20_replicate_study["CI_ub_BCF_SSL", ] > BCF_FULL_DATA[2])))
IPW_BART_CC_cover = c(cover_80 = mean((missing_80_replicate_study["CI_lb_IPW_BART_CC", ] < BCF_FULL_DATA[2]) & 
                                         (missing_80_replicate_study["CI_ub_IPW_BART_CC", ] > BCF_FULL_DATA[2])), 
                       cover_50 = mean((missing_50_replicate_study["CI_lb_IPW_BART_CC", ] < BCF_FULL_DATA[2]) & 
                                         (missing_50_replicate_study["CI_ub_IPW_BART_CC", ] > BCF_FULL_DATA[2])), 
                       cover_20 = mean((missing_20_replicate_study["CI_lb_IPW_BART_CC", ] < BCF_FULL_DATA[2]) & 
                                         (missing_20_replicate_study["CI_ub_IPW_BART_CC", ] > BCF_FULL_DATA[2])))
IPW_BART_SSL_cover = c(cover_80 = mean((missing_80_replicate_study["CI_lb_IPW_BART_SSL", ] < BCF_FULL_DATA[2]) & 
                                   (missing_80_replicate_study["CI_ub_IPW_BART_SSL", ] > BCF_FULL_DATA[2])), 
                 cover_50 = mean((missing_50_replicate_study["CI_lb_IPW_BART_SSL", ] < BCF_FULL_DATA[2]) & 
                                   (missing_50_replicate_study["CI_ub_IPW_BART_SSL", ] > BCF_FULL_DATA[2])), 
                 cover_20 = mean((missing_20_replicate_study["CI_lb_IPW_BART_SSL", ] < BCF_FULL_DATA[2]) & 
                                   (missing_20_replicate_study["CI_ub_IPW_BART_SSL", ] > BCF_FULL_DATA[2])))

# export results to CSV for safekeeping
nhefs_sims_80_file = file.path(project_dir, "outputs", "nhefs_sims_80.csv")
nhefs_sims_50_file = file.path(project_dir, "outputs", "nhefs_sims_50.csv")
nhefs_sims_20_file = file.path(project_dir, "outputs", "nhefs_sims_20.csv")
write.csv(missing_80_replicate_study, nhefs_sims_80_file)
write.csv(missing_50_replicate_study, nhefs_sims_50_file)
write.csv(missing_20_replicate_study, nhefs_sims_20_file)

# close log file
close(output_file)