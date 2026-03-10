# --- 1. Dependencies ---
library(dplyr)
library(fastDummies)
library(glmnet)
library(lmtest)
library(sandwich)

# --- 2. The PDS Engine Function ---
run_pds <- function(df, y_name, d_name, x_vars, family = "gaussian") {
  
  Y <- df[[y_name]]
  D <- ifelse(df[[d_name]] == "Advisory", 1, 0)
  X_mat <- as.matrix(df[, x_vars])
  
  if (family == "gaussian") {
    # --- Standard Linear PDS ---
    cv_d <- cv.glmnet(X_mat, D, family = "gaussian", standardize = TRUE, nfolds = 10)
    vars_d <- rownames(coef(cv_d, s = "lambda.min"))[which(as.numeric(coef(cv_d, s = "lambda.min")) != 0)][-1]
    
    cv_y <- cv.glmnet(X_mat, Y, family = "gaussian", standardize = TRUE, nfolds = 10)
    vars_y <- rownames(coef(cv_y, s = "lambda.min"))[which(as.numeric(coef(cv_y, s = "lambda.min")) != 0)][-1]
    
  } else {
    # --- Weighted Logit PDS (Your Procedure) ---
    # 1. Selection for Y (Logit LASSO)
    cv_y <- cv.glmnet(X_mat, Y, family = "binomial", type.measure = "deviance", nfolds = 10)
    vars_y <- rownames(coef(cv_y, s = "lambda.min"))[which(as.numeric(coef(cv_y, s = "lambda.min")) != 0)][-1]
    
    # 2. Generate Weights for the D-selection
    # We fit a quick logit to get the density weights (w)
    tmp_df <- data.frame(Y = Y, X_mat[, vars_y, drop = FALSE])
    fit_y <- glm(Y ~ ., data = tmp_df, family = "binomial")
    w <- dlogis(predict(fit_y, type = "link"))
    
    # 3. Selection for D (Weighted Linear LASSO)
    cv_d <- cv.glmnet(X_mat, D, weights = w, family = "gaussian", nfolds = 10)
    vars_d <- rownames(coef(cv_d, s = "lambda.min"))[which(as.numeric(coef(cv_d, s = "lambda.min")) != 0)][-1]
  }
  
  # Step 3: Union and Post-Lasso
  selected_vars <- unique(c(vars_d, vars_y))
  selected_vars <- setdiff(selected_vars, c("(Intercept)", y_name, d_name))
  
  final_df <- data.frame(Y = Y, D = D, df[, selected_vars, drop = FALSE])
  post_model <- glm(Y ~ ., data = final_df, family = family)
  
  robust_results <- coeftest(post_model, vcov = vcovHC(post_model, type = "HC3"))
  return(list(results = robust_results, count = length(selected_vars)))
}

# --- 3. Pre-processing ---
# Same cleaning logic as your IPW/ATO code
mydf <- read.csv("/Synthetic_Replication_Data.csv", stringsAsFactors = TRUE)
mydf[which(mydf$YearStart == 2015),'YearStart'] <- 2016
mydf$MonthStart <- factor(mydf$MonthStart)
mydf$YearStart <- factor(mydf$YearStart)
mydf<- mydf[ ,  -which(names(mydf) %in% c('MonthYearStart'))]

# Dummy encoding factors
fact_cols <- names(Filter(is.factor, mydf))
fact_cols <- fact_cols[fact_cols != "RequestedServices"]
mydf <- dummy_cols(mydf, select_columns = fact_cols, remove_first_dummy = TRUE, remove_selected_columns = TRUE)

# --- 4. Loop through Outcomes ---
continuous_outcomes <- c('InvestmentValue', 'DebtEnd')
binary_outcomes <- c('ImpactResults', 'ImpactInvestment', 'IfBadArrearsMove', 'IfBankedMove', 'IfGoodArrearsMove')

for(y in continuous_outcomes) {
  print(paste("### PDS Results for:", y))
  print(run_pds(mydf, y, "RequestedServices", setdiff(names(mydf), c("RequestedServices", y_vars)), family = "gaussian"))
}

for(y in binary_outcomes) {
  print(paste("### PDS Results for:", y))
  print(run_pds(mydf, y, "RequestedServices", setdiff(names(mydf), c("RequestedServices", y_vars)), family = "binomial"))
}

# --- 5. Falsification Test (Unconfoundedness Check) ---
# We use DebtStart as the outcome. 
# It should NOT be predicted by RequestedServices if the model is balanced.

y_falsification <- "DebtStart"

# Ensure DebtStart is NOT in your covariate list (x_vars) 
# and that other future outcomes are also excluded to prevent leakage.
x_vars_falsification <- setdiff(names(mydf), c("RequestedServices", y_vars, y_falsification))

print(paste("### PDS Falsification Test for:", y_falsification))

# Run PDS with Gaussian family since DebtStart is continuous
falsification_results <- run_pds(
  df = mydf, 
  y_name = y_falsification, 
  d_name = "RequestedServices", 
  x_vars = x_vars_falsification, 
  family = "gaussian"
)

print(falsification_results$results)

