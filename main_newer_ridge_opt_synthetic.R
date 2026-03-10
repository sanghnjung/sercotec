# --- 1. Dependencies and Data Loading ---
library(dplyr)
library(fastDummies)
library(glmnet)
library(doParallel)
library(PSweight)

# Load synthetic data
mydf <- read.csv("C:/Users/jshsa/Desktop/project/project1/sercotec/Synthetic_Replication_Data.csv", stringsAsFactors = TRUE)
mydf[which(mydf$YearStart == 2015),'YearStart'] <- 2016
mydf$MonthStart <- factor(mydf$MonthStart)
mydf$YearStart <- factor(mydf$YearStart)
mydf<- mydf[ ,  -which(names(mydf) %in% c('MonthYearStart'))]
# --- 2. Helper Functions ---

# Ridge Propensity Score Fitting
fit_ps_ridge <- function(df, xvars, zname, trtgrp, s = "lambda.min") {
  d <- as.integer(df[[zname]] == trtgrp)
  Xmat <- as.matrix(df[, xvars, drop = FALSE])
  
  cv <- cv.glmnet(Xmat, d, family = "binomial", alpha = 0, standardize = TRUE)
  ps <- as.numeric(predict(cv, s = s, newx = Xmat, type = "response"))
  ps <- pmax(pmin(ps, 0.9999), 0.0001)
  return(ps)
}

# Optimal Trimming Threshold (Crump et al., 2009)
find_opt_delta <- function(ps, grid = seq(0.01, 0.25, by = 0.005)) {
  obj <- sapply(grid, function(d) {
    keep <- (ps > d) & (ps < 1 - d)
    pk <- mean(keep)
    if (pk == 0) return(Inf)
    mean((1 / (ps * (1 - ps))) * keep, na.rm = TRUE) / (pk^2)
  })
  grid[which.min(obj)]
}

# Calculate Causal Estimates (ATE, ATT, ATC, ATO)
get_estimates <- function(df, ps, trtgrp, yname, zname, opt_delta) {
  # Trimmed Sample for IPW
  keep <- ps > opt_delta & ps < 1 - opt_delta
  df_t <- df[keep, ]
  ps_t <- ps[keep]
  
  n_t <- length(ps_t)
  
  trt_t <- (df_t[[zname]] == trtgrp)
  
  # Full Sample for Overlap Weights (ATO)
  trt_f <- (df[[zname]] == trtgrp)
  ps_f  <- ps
  
  results <- list()
  
  for (y in yname) {
    Y_t <- df_t[[y]]
    Y_f <- df[[y]]
    
    # Normalized IPW Weights
    w1_ate <- 1 / ps_t; w0_ate <- 1 / (1 - ps_t)
    w1_att <- w1_att <- rep(1, n_t);        w0_att <- ps_t / (1 - ps_t)
    w1_atc <- (1 - ps_t) / ps_t; w0_atc <- rep(1, n_t)
    
    # Overlap Weights
    w1_ato <- 1 - ps_f; w0_ato <- ps_f
    
    calc_norm <- function(Y, trt, w1, w0) {
      (sum(Y[trt] * w1[trt]) / sum(w1[trt])) - (sum(Y[!trt] * w0[!trt]) / sum(w0[!trt]))
    }
    
    results[[y]] <- c(
      ATE = calc_norm(Y_t, trt_t, w1_ate, w0_ate),
      ATT = calc_norm(Y_t, trt_t, w1_att, w0_att),
      ATC = calc_norm(Y_t, trt_t, w1_atc, w0_atc),
      ATO = calc_norm(Y_f, trt_f, w1_ato, w0_ato)
    )
  }
  return(unlist(results))
}

# --- 3. Pre-processing ---
# Define variables
y_vars <- c('ImpactInvestment', 'InvestmentValue', 'ImpactResults', 'DebtEnd', 
            'IfBadArrearsMove', 'IfBankedMove', 'IfGoodArrearsMove')
z_name <- 'RequestedServices'
trt_grp <- 'Advisory'

# Clean and Dummy Encode
X <- mydf
X$RequestedServices <- factor(X$RequestedServices, levels = c('Training', 'Advisory'))

# Identify factor columns for dummying (excluding treatment)
fact_cols <- names(X)[sapply(X, is.factor)]
fact_cols <- fact_cols[fact_cols != z_name]
X <- dummy_cols(X, select_columns = fact_cols, remove_first_dummy = TRUE, remove_selected_columns = TRUE)

x_vars_all <- setdiff(names(X), c(z_name, y_vars))

# --- 4. Main Analysis ---
ps_main <- fit_ps_ridge(X, x_vars_all, z_name, trt_grp)
opt_delta <- find_opt_delta(ps_main)
main_est <- get_estimates(X, ps_main, trt_grp, y_vars, z_name, opt_delta)

# --- 5. Bootstrap Parallelization ---
B <- 10
nodes <- unlist(strsplit(Sys.getenv("NODESLIST"), " "))
cl <- makeCluster(length(nodes), type = "PSOCK")
registerDoParallel(cl)

# Export essential items to workers
clusterExport(cl, c("X", "x_vars_all", "trt_grp", "y_vars", "z_name", "opt_delta", 
                    "fit_ps_ridge", "get_estimates", "find_opt_delta"))
clusterEvalQ(cl, { library(glmnet); library(dplyr) })
clusterSetRNGStream(cl, 6262)

boot_results <- foreach(i = 1:B, .combine = cbind, .errorhandling = 'remove') %dopar% {
  # Stratified Bootstrap
  idx <- c(sample(which(X[[z_name]] == trt_grp), replace = TRUE),
           sample(which(X[[z_name]] != trt_grp), replace = TRUE))
  df_b <- X[idx, ]
  
  ps_b <- fit_ps_ridge(df_b, x_vars_all, z_name, trt_grp)
  get_estimates(df_b, ps_b, trt_grp, y_vars, z_name, opt_delta)
}

stopCluster(cl)

# --- 6. Results Formatting ---
calc_stats <- function(boot_mat, point_est) {
  # Calculate Standard Deviation from the bootstrap distribution
  sds <- apply(boot_mat, 1, sd, na.rm = TRUE)
  
  # Calculate 95% Confidence Intervals using the Percentile Method
  ci_lb <- apply(boot_mat, 1, quantile, probs = 0.025, na.rm = TRUE)
  ci_rb <- apply(boot_mat, 1, quantile, probs = 0.975, na.rm = TRUE)
  
  # Removed p-value calculation as requested
  data.frame(
    Estimate = point_est, 
    StdErr = sds, 
    LCI = ci_lb, 
    UCI = ci_rb
  )
}

final_tab <- calc_stats(boot_results, main_est)
print(final_tab)

# --- 7.1. falsification Pre-processing ---
y_placebo <- 'DebtStart'

# Create a fresh dataset for placebo (Removing future outcomes to prevent leakage)
X_placebo <- mydf
X_placebo$RequestedServices <- factor(X_placebo$RequestedServices, levels = c('Training', 'Advisory'))

# Identify and Dummy factors
fact_cols_p <- names(X_placebo)[sapply(X_placebo, is.factor)]
fact_cols_p <- fact_cols_p[fact_cols_p != z_name]
X_placebo <- dummy_cols(X_placebo, select_columns = fact_cols_p, 
                        remove_first_dummy = TRUE, remove_selected_columns = TRUE)

# Filter out the OTHER outcome variables so the model only sees DebtStart and Covariates
# We keep only DebtStart as the 'y' and everything else as 'x'
x_vars_placebo <- setdiff(names(X_placebo), c(z_name, y_vars, y_placebo))

# --- 7.2. Main Placebo Analysis ---
ps_placebo <- fit_ps_ridge(X_placebo, x_vars_placebo, z_name, trt_grp)
opt_delta_p <- find_opt_delta(ps_placebo)
placebo_est <- get_estimates(X_placebo, ps_placebo, trt_grp, y_placebo, z_name, opt_delta_p)

# --- 7.3. Placebo Bootstrap Parallelization ---
cl_p <- makeCluster(length(nodes), type = "PSOCK")
registerDoParallel(cl_p)

clusterExport(cl_p, c("X_placebo", "x_vars_placebo", "trt_grp", "y_placebo", "z_name", "opt_delta_p", 
                      "fit_ps_ridge", "get_estimates", "find_opt_delta"))
clusterEvalQ(cl_p, { library(glmnet); library(dplyr) })
clusterSetRNGStream(cl_p, 6262)

boot_placebo <- foreach(i = 1:B, .combine = cbind, .errorhandling = 'remove') %dopar% {
  idx <- c(sample(which(X_placebo[[z_name]] == trt_grp), replace = TRUE),
           sample(which(X_placebo[[z_name]] != trt_grp), replace = TRUE))
  df_b <- X_placebo[idx, ]
  
  ps_b <- fit_ps_ridge(df_b, x_vars_placebo, z_name, trt_grp)
  get_estimates(df_b, ps_b, trt_grp, y_placebo, z_name, opt_delta_p)
}

stopCluster(cl_p)

# --- 7.4. Placebo Results Formatting ---
placebo_tab <- calc_stats(boot_placebo, placebo_est)
print("Placebo Check Results (Should be non-significant):")
print(placebo_tab)