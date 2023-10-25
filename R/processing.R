

#' Extract DFC with Adaptive Lasso
#'
#' @param data scaled count gene x cell matrix
#' @param label logical or binary vector labeling target clusters
#' @param gamma Parameter to control the effect of penalty
#' @param lambda_penalty Parameter to control which lambda is used to calculate
#' penalty
#' @param lambda_weight Parameter to control which lambda is used to calculate
#' weights
#' @param seed seed
#'
#' @returns \item{Ridge}{The model of ridge regression to calculate penalry.}
#' \item{AdaLasso}{The model of Adaptive Lasso for DFC extraction.}
#' \item{weights}{A data frame of extracted features and the weights. When
#' return_Model set FALSE, only weigths are returned.}
#'
#' @importFrom glmnet cv.glmnet
#' @importFrom glmnet coef.glmnet
#'
AdaLasso <- function(
    data, label, gamma = 1, lambda_penalty = "1se", lambda_weight = "1se",
    seed = NULL
    ) {
  lambda_penalty <- switch(
    lambda_penalty,
    "1se" = "lambda.1se",
    "min" = "lambda.min"
  )
  lambda_weight <- switch(
    lambda_weight,
    "1se" = "lambda.1se",
    "min" = "lambda.min"
  )
  if(is.null(seed)) seed <- 0
  mat <- as.matrix(scale(t(data)))
  #mat <- Matrix::drop0(mat)
  cat("Calculating penalties...\n")
  set.seed(seed)
  ridge <-  cv.glmnet(
    mat, label,
    family = "binomial", alpha = 0,
    lambda = exp(seq(-7,1,len=40))
  )
  ridge_coef <- coef.glmnet(ridge, s = lambda_penalty)[-1]
  penalty <- 1/abs(ridge_coef)^gamma
  cat("Calculating weights...\n")
  alasso <- cv.glmnet(
    mat, label,
    family = "binomial", alpha = 1,
    penalty = penalty
  )
  adalass_coef <- coef.glmnet(alasso, s = lambda_weight)[-1]
  res <- data.frame(
    feature = colnames(mat),
    weight = adalass_coef
  )
  res <- res[res$weight!=0,]
  cat(paste0(nrow(res)," features were selected.\n"))
  models <- list(
    Ridge = ridge, AdaLasso = alasso,
    weights = res
    )
  class(models) <- 'dfc_models'
  return(models)
}

#' Sure independence screening for pre-screening of features.
#'
#' @param data scaled count gene x cell matrix
#' @param label logical or binary vector labeling target clusters
#' @param max_feature maximum number to extract
#' @param min_feature minimum number to extract
#'
#' @returns screened data
#'
sis <- function(
    data, label, min_feature = NULL, max_feature = NULL
    ) {
  nsis <- floor(ncol(data)/(4*log(ncol(data))))
  if(!is.null(min_feature)) nsis <- max(min_feature,nsis)
  if(!is.null(max_feature)) nsis <- min(max_feature,nsis)
  coefs <- abs(cor(t(data), as.integer(label)))
  ordercoef = order(coefs, decreasing = TRUE)
  sis_data <- data[ordercoef[1:nsis],]
  cat(paste0(nsis," features were selected by SIS.\n"))
  return(sis_data)
}
