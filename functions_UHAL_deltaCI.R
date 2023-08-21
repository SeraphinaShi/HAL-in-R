###############################################################################
#'  fit undersoomthed HAL function with global criterion
#'
#' @details The procedure select the lambda that satisfies the global 
#' undersmoothing criterion (\eqn{ P_n(\phi_{s,i}(Y-\bar{Q}_{n,\lambda}))\leq \freq{\sigma_n}{\sqrt{n}log(n)} }).
#' It performs the following steps:
#'     1). get all the basis functions from cv-hal 
#'     2). calculate a grid of new lambdas by scaling the cv-lambda from
#'          cv-hal with seq(from=0, to=-3, length=Nlam)
#'     3). refit lasso using \code{\link[glmnet]{glmnet}}) with all the basis 
#'         functions and with the grid of new lambdas
#'     4).calculate the scores (\eqn{ P_n(\phi_{s,i}(Y-\bar{Q}_{n,\lambda})) })of each lasso fit     
#'     5).find the biggest lambda that max score (\eqn{ \leq \freq{\sigma_n}{\sqrt{n}log(n)} })
#' @param X An input \code{matrix} with dimensions number of observations -by-
#'  number of covariates that will be used to derive the design matrix of basis
#'  functions.
#' @param Y A \code{numeric} vector of the outcome variable of observations. 
#' @param fit_init The initial HAL fit object from the output list of \code{undersmooth_init}.
#' @param Nlam Number of lambda candidates. The sequence ranges from \code{fit_init$lambda_star} to
#' \code{fit_init$lambda_star*10^(-3)}.
#' @param family A \code{character} or a \code{\link[stats]{family}} object
#'  (supported by \code{\link[glmnet]{glmnet}}) specifying the error/link
#'  family for a generalized linear model. See more details in \code{fit_hal}.

undersmooth_hal <- function(X,
                            Y,
                            fit_init,
                            Nlam = 20,
                            family = c("gaussian", "binomial", "poisson", "cox")
                            ){
  
  n = length(Y)
  
  nonzero_col_init = which(fit_init$coefs[-1] != 0)
  if (length(nonzero_col_init) == 0){
    res <- list("lambda_init" = fit_init$lambda_star,
                "lambda_under" = fit_init$lambda_star)
    return(res)
  }
  
  #  refit on a grid of lambda (scaled on the lambda from cv-hal)
  us_lambda <- fit_init$lambda_star*10^seq(from=0, to=-3, length=Nlam)
  us_fit <- glmnet(fit_init$x_basis, Y, lambda=us_lambda, family = family, standardize = FALSE)
  
  if(identical(us_fit$df, 0)){
    res <- list("lambda_init" = fit_init$lambda_star,
                "lambda_under" = fit_init$lambda_star)
    return(res)
  }
  
  
  if (family == "binomial"){
    preds_init <- predict(fit_init, new_data = X, type = "response")
    pred_mat <- predict(us_fit, fit_init$x_basis, type = "response")
  }else {
    preds_init <- predict(fit_init, new_data = X)
    pred_mat <- predict(us_fit, fit_init$x_basis)
  }
  resid_mat <- pred_mat - Y
  
  ## estimates of sd in each direction using initial fit
  resid_init <- preds_init - Y
  
  basis_mat_init <- as.matrix(fit_init$x_basis)
  basis_mat_init <- as.matrix(basis_mat_init[, nonzero_col_init])
  
  sd_est  <- apply(basis_mat_init, 2, function(phi) sd(resid_init*phi))
  
  
  ## calculate scores
  max_score <- get_maxscore(basis_mat = basis_mat_init,
                            resid_mat = resid_mat,
                            sd_est = sd_est,
                            Nlam = Nlam, us_fit = us_fit)
  
  ## get the first lambda that satisfies the criteria
  lambda_under <- us_lambda[max_score <= 1/(sqrt(n)*log(n))][1] # over under-smoothing 
  
  if (is.na(lambda_under)){
    res <- list("lambda_init" = fit_init$lambda_star,
                "lambda_under" = fit_init$lambda_star)
    return(res)
  }
  
  # collect results
  coef_mat <- as.matrix(us_fit$beta)
  
  spec_under <- list("lambda" = us_lambda,
                     "l1_norm" = NA,
                     "n_coef" = NA)
  
  spec_under$l1_norm <- apply(coef_mat, 2, function(x){sum(abs(x))})
  spec_under$n_coef <- apply(coef_mat, 2, function(x){sum(x != 0)})
  
  res <- list("lambda_init" = fit_init$lambda_star,
              "lambda_under" = lambda_under,
              "spec_under" = spec_under)
  return(res)
}


###############################################################################
#'  undersoomthed HAL helper function for global criterion
#'
#' @details For each candidate lambda, do:
#'     1). standardize the score formed by each basis.
#'     2). calculate the mean of the standardized scores for each basis.
#' Select the max of the mean.
#' @param basis_mat The selected basis matrix from initial fit for undersmoothing,
#'  obtained from the output list of \code{undersmooth_init}.
#' @param resid_mat The residual \code{numeric} with each column the residuals correspongding to a lambda.
#' @param sd_est A \code{numeric} vector containing the sd of each column of \code{basis_mat}.
#' @param Nlam Number of lambda candidates.
#' @param us_fit The \code{glmnet} fit of the sequence of candidate lambdas.

get_maxscore <- function(basis_mat, resid_mat, sd_est, Nlam, us_fit){
  
  basis_mat_sd <- sweep(basis_mat, 2, sd_est, FUN = '/')
  score_all <- apply(basis_mat_sd, 2, function(u) {
    score_mat <- resid_mat * u
    score_mean <- apply(score_mat, 2, mean)
  })
  # absolute value
  max_score <- apply(abs(score_all), 1, max, na.rm=T)
  return(max_score)
}



###############################################################################
#'  calculate ATE estimate and a confidence interval using delta method.
#'
#' @details This is an inference function for average treatment effect (ATE). 
#'  With provided X, Y, and their corresponding HAL fit, it will calculate the 
#'  ATE estimate, SE, and confidence intervals using delta method. 
#' @param X An input \code{matrix} with dimensions number of observations -by-
#'  number of covariates that will be used to derive the design matrix of basis
#'  functions.The same \code{X} used in \code{fit_hal} function.
#' @param Y A \code{numeric} vector of the outcome variable of observations. 
#'  The same \code{Y} used in \code{fit_hal} function.
#' @param treatment A \code{character} specifying the treatment column name in \code{X}
#' @param condition A \code{character} specifying the condition column name in 
#'  \code{X}, only used when estimating the CATE.
#' @param condition_level A \code{character} or a \code{numeric} specifying the 
#'  condition level, only used when estimating the CATE.
#' @param hal_fit The HAL fit object from the output of \code{fit_hal} function.
#' @param alpha The statistical significance level for the influence curve.  
#' @param family A \code{character} or a \code{\link[stats]{family}} object
#'  (supported by \code{\link[glmnet]{glmnet}}) specifying the error/link
#'  family for a generalized linear model. See more details in \code{fit_hal}.
ATE_inferece <- function(X, Y, treatment, condition = NA, condition_level = NA, hal_fit, alpha = 0.05, family = "gaussian"){
  
  n = length(Y)
  treat_idx = which(colnames(X) == treatment)
  
  X = as.matrix(X)
  
  if(! is.na(condition)){
    cond_idx = which(colnames(X) == condition)
    obj_cond_idx = X[, condition] == condition_level
    X_c = X[obj_cond_idx, ]
    Y_c = Y[obj_cond_idx]
  } else {
    X_c = X
    Y_c = Y
  }
  
  # 1. get estimate
  ## intervene to set A=a and generate the counterfactual outcomes Y.a
  X.1 = as.matrix(X_c) 
  X.1[,treat_idx] = 1
  Y.1 = predict(hal_fit, new_data = X.1)
  
  X.0 = as.matrix(X_c) 
  X.0[,treat_idx] = 0
  Y.0 = predict(hal_fit, new_data = X.0)
  
  ## Estimation 
  ATE_est = mean(Y.1 - Y.0)
  
  # 2. get inferences
  ## 2.1. influence curve of beta
  coef <- hal_fit$coefs
  basis_mat <- cbind(1, as.matrix(hal_fit$x_basis))
  
  nonzero_idx <- which(coef != 0)
  
  coef_nonzero <- coef[nonzero_idx]
  basis_mat_nonzero <- as.matrix(basis_mat[, nonzero_idx])
  
  if(family == 'binomial'){
    Y_hat = predict(hal_fit, new_data = X, type = "response")
  } else {
    Y_hat = predict(hal_fit, new_data = X)
  }
  
  IC_beta = cal_IC_for_beta(X = basis_mat_nonzero,
                            Y = Y,
                            Y_hat = Y_hat,
                            beta_n = coef_nonzero,
                            family = family)
  
  if(! is.na(condition)){
    IC_beta = IC_beta[, obj_cond_idx]
  }
  
  ## 2.2. influence curve of point estimate
  x_basis.1 <- make_design_matrix(X.1, hal_fit$basis_list, p_reserve = 0.75)
  x_basis_nonzero.1 <- as.matrix(cbind(1, x_basis.1)[, nonzero_idx])
  
  IC_EY.1 <- cal_IC_for_EY(X = x_basis_nonzero.1, 
                           beta_n = coef_nonzero, 
                           IC_beta = IC_beta,
                           family = family)
  
  x_basis.0 <- make_design_matrix(X.0, hal_fit$basis_list, p_reserve = 0.75)
  x_basis_nonzero.0 <- as.matrix(cbind(1, x_basis.0)[, nonzero_idx])
  
  IC_EY.0 <- cal_IC_for_EY(X = x_basis_nonzero.0, 
                           beta_n = coef_nonzero, 
                           IC_beta = IC_beta,
                           family = family)
  
  IC_ATE = IC_EY.1 - IC_EY.0
  
  ## 2.3. standard error & confidence intervals
  SE_ATE = sqrt(var(IC_ATE)/n)
  
  ci_lwr_ATE <- ATE_est - abs(qnorm(alpha/2)) * SE_ATE
  ci_upr_ATE <- ATE_est + abs(qnorm(alpha/2)) * SE_ATE
  
  if(is.na(condition)){
    Target_param = paste0("E[Y_1] - E[Y_0]")
    returns = list(Target_parameter = paste0("E[Y_1] - E[Y_0]"),
                   Estimate = ATE_est,
                   SE = SE_ATE,
                   CI = c(ci_lwr_ATE, ci_upr_ATE))
  } else {
    Target_param = paste0("E[Y_1|", condition, "=", condition_level, "] - E[Y_0|", condition, "=", condition_level, "]")
  }
  
  returns = list(Target_parameter = Target_param,
                 Estimate = ATE_est,
                 SE = SE_ATE,
                 CI = c(ci_lwr_ATE, ci_upr_ATE))
  
  return(returns)

}

###############################################################################
#'  calculate the counterfactual outcome mean estimate and a confidence interval using delta method.
#'
#' @details This is an inference function for counterfactual outcome mean. 
#'  With provided X, Y, and their corresponding HAL fit, it will calculate the 
#'  counterfactual outcome mean estimate, SE, and confidence intervals using delta method. 
#' @param X An input \code{matrix} with dimensions number of observations -by-
#'  number of covariates that will be used to derive the design matrix of basis
#'  functions.The same \code{X} used in \code{fit_hal} function.
#' @param Y A \code{numeric} vector of the outcome variable of observations. 
#'  The same \code{Y} used in \code{fit_hal} function.
#' @param treatment A \code{character} specifying the treatment column name in \code{X}
#' @param treatment_level A \code{character} or a \code{numeric} value specifying 
#'  the counterfactual treatment level.
#' @param condition A \code{character} specifying the condition column name in 
#'  \code{X}, only used when estimating the CATE.
#' @param condition_level A \code{character} or a \code{numeric} specifying the 
#'  condition level, only used when estimating the CATE.
#' @param hal_fit The HAL fit object from the output of \code{fit_hal} function.
#' @param alpha The statistical significance level for the influence curve.  
#' @param family A \code{character} or a \code{\link[stats]{family}} object
#'  (supported by \code{\link[glmnet]{glmnet}}) specifying the error/link
#'  family for a generalized linear model. See more details in \code{fit_hal}.
CounterfactualMean_inferece <- function(X, Y, treatment, treatment_level, condition = NA, condition_level = NA, hal_fit, alpha = 0.05, family = "gaussian"){
  
  n = length(Y)
  treat_idx = which(colnames(X) == treatment)
  
  if(! is.na(condition)){
    cond_idx = which(colnames(X) == condition)
    obj_cond_idx = X[, condition] == condition_level
    X_c = X[obj_cond_idx, ]
    Y_c = Y[obj_cond_idx]
  } else {
    X_c = X
    Y_c = Y
  }
  
  # 1. get estimate
  ## intervene to set A=a and generate the counterfactual outcomes Y.a
  X.a = as.matrix(X_c) 
  X.a[,treat_idx] = treatment_level
  Y.a = predict(hal_fit, new_data = X.a)
  
  ## Estimation 
  CounterfactualMean_est = mean(Y.a)
  
  # 2. get inferences
  ## 2.1. influence curve of beta
  coef <- hal_fit$coefs
  basis_mat <- cbind(1, as.matrix(hal_fit$x_basis))
  
  nonzero_idx <- which(coef != 0)
  
  coef_nonzero <- coef[nonzero_idx]
  basis_mat_nonzero <- as.matrix(basis_mat[, nonzero_idx])
  
  if(family == 'binomial'){
    Y_hat = predict(hal_fit, new_data = X, type = "response")
  } else {
    Y_hat = predict(hal_fit, new_data = X)
  }
  
  IC_beta = cal_IC_for_beta(X = basis_mat_nonzero,
                            Y = Y,
                            Y_hat = Y_hat,
                            beta_n = coef_nonzero,
                            family = family)
  
  if(! is.na(condition)){
    IC_beta = IC_beta[, obj_cond_idx]
  }
  
  ## 2.2. influence curve of point estimate
  x_basis.a <- make_design_matrix(X.a, hal_fit$basis_list, p_reserve = 0.75)
  x_basis_nonzero.a <- as.matrix(cbind(1, x_basis.a)[, nonzero_idx])
  
  IC_EY.a <- cal_IC_for_EY(X = x_basis_nonzero.a, 
                           beta_n = coef_nonzero, 
                           IC_beta = IC_beta,
                           family = family)
  
  ## 2.3. standard error & confidence intervals
  SE = sqrt(var(IC_EY.a)/n)
  
  ci_lwr <- CounterfactualMean_est - abs(qnorm(alpha/2)) * SE
  ci_upr <- CounterfactualMean_est + abs(qnorm(alpha/2)) * SE
  

  
  if(is.na(condition)){
    returns = list(Target_parameter = paste0("EE[Y|A=", treatment_level, ", W]"),
                   Estimate = CounterfactualMean_est,
                   SE = SE,
                   CI = c(ci_lwr, ci_upr))
  } else {
    returns = list(Target_parameter = paste0("EE[Y|A=", treatment_level, ",", ondition, "=", condition_level, " W]"),
                   Estimate = CounterfactualMean_est,
                   SE = SE,
                   CI = c(ci_lwr, ci_upr))
  }
  
  return(returns)
  
}

###############################################################################
#'  HAL based inferences helper function for calculating the influence curve of beta
#'
#' @details Calculate the influence curve of beta, $D_p(\Phi,Y)$
#' 
#' @param X The basis matrix with non-zero coefficients from initial fit.
#' @param Y A \code{numeric} vector of the outcome variable of observations. 
#'  The same \code{Y} used in \code{fit_hal} function.
#' @param Y_hat A \code{numeric} vector of the predicted outcome variable of 
#'  observations from \code{predict} function using on the HAL fit object. 
#' @param beta_n A \code{numeric} vector containing the non-zero coefficients.
#' @param family A \code{character} or a \code{\link[stats]{family}} object
#'  (supported by \code{\link[glmnet]{glmnet}}) specifying the error/link
#'  family for a generalized linear model. See more details in \code{fit_hal}.

cal_IC_for_beta <- function(X, Y, Y_hat, beta_n, family = 'binomial'){
  n <- dim(X)[1] 
  p <- length(beta_n)
  if (!is.matrix(X)) X <- as.matrix(X)
  
  # 1. calculate score: X'(Y - phi(X))
  res <- Y-Y_hat
  score <- sweep(t(X), 2, res, `*`)
  
  # 2. calculate d_phi, the derivative of phi:
  if(family == 'binomial'){
    d_phi_scaler <- as.vector(exp(- beta_n %*% t(X)) / ((1 + exp(- beta_n %*% t(X)))^2))
    d_phi <- sweep(X, 1, d_phi_scaler, `*`)
  } else {
    d_phi = - X
  }
  
  # 3. -E_{P_n}(X d_phi)^(-1)
  tmat <- t(X) %*% d_phi / n
  if(! is.matrix(try(solve(tmat), silent = TRUE))){
    return(NA)
  }
  tmat <- -solve(tmat)
  
  # 4. calculate influence curves
  IC <- tmat %*% score
  
  return(IC)
}


###############################################################################
#'  HAL based inferences helper function for calculating the influence curve of Y
#'
#' @details Calculate the influence curve of Y, 
#'  $IC_{\psi} = \frac{d}{d\beta}\psi_{P_n}(\Phi) * D_p(\Phi,Y) $
#' 
#' @param X A basis matrix with non-zero coefficients from initial fit.
#' @param beta_n A \code{numeric} vector containing the non-zero coefficients.
#' @param IC_beta A \code{numeric} vector containing the influence curves of beta.
#' @param family A \code{character} or a \code{\link[stats]{family}} object
#'  (supported by \code{\link[glmnet]{glmnet}}) specifying the error/link
#'  family for a generalized linear model. See more details in \code{fit_hal}.
cal_IC_for_EY <- function(X, beta_n, IC_beta, family = 'binomial'){
  if (!is.matrix(X)) X <- as.matrix(X)
  
  if(family == 'binomial'){
    d_phi_scaler_new <- as.vector(exp(- beta_n %*% t(X)) / ((1 + exp(- beta_n %*% t(X)))^2))
    d_phi_new <- sweep(X, 1, d_phi_scaler_new, `*`)
  } else {
    d_phi_new = X
  }
  
  IC = diag(d_phi_new %*% IC_beta)
  
  return(IC)
}








