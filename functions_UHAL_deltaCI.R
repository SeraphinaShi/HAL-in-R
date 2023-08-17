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
#' @param Y A \code{numeric} vector of observations of the outcome variable.
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
  
  
  if (family != "binomial"){
    preds_init <- predict(fit_init, new_data = X)
    pred_mat <- predict(us_fit, fit_init$x_basis)
  }else {
    preds_init <- predict(fit_init, new_data = X, type = "response")
    pred_mat <- predict(us_fit, fit_init$x_basis, type = "response")
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
#' @param resid_mat The residual matrix with each column the residuals correspongding to a lambda.
#' @param sd_est A numeric vector containing the sd of each column of \code{basis_mat}.
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