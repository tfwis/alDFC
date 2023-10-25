#' Extract DFC with Adaptive Lasso
#'
#' @param data scaled count gene x cell matrix
#' @param label logical or binary vector labeling target clusters
#' @param gamma Parameter to control the effect of penalty
#' @param return_Model Return regression models or only weights
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
    data, label, gamma = 1,
    return_Model = TRUE,
    lambda_penalty = "1se",
    lambda_weight = "1se",
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
  if (return_Model) {
    models <- list(
      Ridge = ridge, AdaLasso = alasso,
      weights = res
      )
    class(models) <- 'dfc_models'
    return(models)
  }else{
    return(res)
  }
}

#' Sure independence screening for pre-screening of features.
#'
#' @param data scaled count gene x cell matrix
#' @param label logical or binary vector labeling target clusters
#' @param min_feature minimum number to extract
#' @param max_feature maximum number to extract
#'
#' @returns screened data
#'
sis <- function(data, label, min_feature = NULL, max_feature = NULL)
{
  nsis <- floor(ncol(data)/(4*log(ncol(data))))
  if(!is.null(min_feature)) nsis <- max(min_feature,nsis)
  if(!is.null(max_feature)) nsis <- min(max_feature,nsis)
  coefs <- abs(cor(t(data), as.integer(label)))
  ordercoef = order(coefs, decreasing = TRUE)
  sis_data <- data[ordercoef[1:nsis],]
  cat(paste0(nsis," features were selected by SIS.\n"))
  return(sis_data)
}

#' Extract DFC
#'
#' @param data SeuratObject or gene x cell Matrix
#' @param ... ...
#' @returns The result of DFC
#'
#' @export
#'
dfc <- function(data, ...) {
  UseMethod(generic = "dfc",object = data)
}

#' Extract DFC
#'
#' @param data SeuratObject
#' @param target_clusters The target cluster number
#' @param assay The assay to use. In default, the return of DefaultAssay() is 
#' used
#' @param gamma The parameter to control the effect of penalty
#' @param return_Model Return regression models or only weights
#' @param seed seed
#' @param lambda_penalty Parameter to control which lambda is used to calculate
#' penalty, "min" or "1se". "min" : binomial deviance is minimized.
#' "1se" : the number of extracted features are minimized within error range of
#' "min".
#' @param lambda_weight Parameter to control which lambda is used to calculate
#' weights.
#' @param SIS Perform screening by SIS or not
#' @param min_feature minimum number to extract
#' @param max_feature maximum number to extract
#'
#' @returns \item{Ridge}{The model of ridge regression to calculate penalry.}
#' \item{AdaLasso}{The model of Adaptive Lasso for DFC extraction.}
#' \item{weights}{A data frame of extracted features and the weights. When
#' return_Model set FALSE, only weigths are returned.}
#'
#' @importFrom Seurat DefaultAssay
#' @importFrom Seurat GetAssay
#' @importFrom Seurat GetAssayData
#' @importFrom Seurat %||%
#'
dfc.Seurat <- function(
    data, target_clusters, assay = NULL,
    gamma = 1, return_Model = FALSE, seed = NULL,
    lambda_penalty = "1se", lambda_weight = "1se",
    SIS = TRUE, min_feature = NULL, max_feature = NULL
    ) {
  cat("Preprocessing...\n")
  assay <- assay %||% DefaultAssay(data)
  smat <- GetAssayData(GetAssay(data,assay),'scale.data')
  if(length(smat)==0) {
    stop("Please, run Seurat::ScaleData")
  }
  target_label <- data$seurat_clusters %in% target_clusters
  if(SIS) {
    smat <- sis(smat, target_label, min_feature, max_feature)
  }
  res <- AdaLasso(
    data = smat, label = target_label,
    gamma = gamma, return_Model = return_Model, seed = seed,
    lambda_penalty = lambda_penalty, lambda_weight = lambda_weight
  )
  return(res)
}

#' Extract DFC
#'
#' @param data scaled count matrix
#' @param target_clusters The target cluster number
#' @param cluster_label Logical or binary vector labeling target clusters
#' @param gamma Parameter to control the effect of penalty
#' @param return_Model Return regression models or only weights
#' @param seed seed
#' @param lambda_penalty Parameter to control which lambda is used to calculate
#' penalty, "min" or "1se". "min" : binomial deviance is minimized.
#' "1se" : the number of extracted features are minimized within error range of
#' "min".
#' @param lambda_weight Parameter to control which lambda is used to calculate
#' weights.
#' @param SIS Perform screening by SIS or not
#' @param min_feature minimum number to extract
#' @param max_feature maximum number to extract
#'
#' @returns \item{Ridge}{The model of ridge regression to calculate penalry.}
#' \item{AdaLasso}{The model of Adaptive Lasso for DFC extraction.}
#' \item{weights}{A data frame of extracted features and the weights. When
#' return_Model set FALSE, only weigths are returned.}
#'
dfc.matrix <- function(
    data, target_clusters, cluster_label,
    gamma = 1, return_Model = FALSE, seed = NULL,
    SIS = TRUE, min_feature = NULL, max_feature = NULL,
    lambda_penalty = "1se", lambda_weight = "1se"
    ) {
  if(ncol(data)!=length(target_label)){
    stop("The column numbers and the label lengths must match.")
  }
  cat("Preprocessing...\n")
  target_label <- cluster_label %in% target_clusters
  if(SIS) {
    data <- sis(data, target_label, min_feature, max_feature)
  }
  res <- AdaLasso(
    data = data, label = target_label,
    gamma = gamma, return_Model = return_Model, seed = seed,
    lambda_penalty = lambda_penalty, lambda_weight = lambda_weight
  )
  return(res)
}

#' Specificity based DFC classification
#'
#' @param data The same one as when `dfc()` is run.#'
#' @param ... ...
#' @returns Thr result data frame with DFC class.
#'
#' @export
#'
dfc_classify <- function(data,...) {
  UseMethod(generic = "dfc_classify",object = data)
}

#' Specificity based DFC classification
#'
#' @param data The same one as when `dfc()` is run.
#' @param dfc_res A result of `dfc()`.
#' @param assay The assay to use. In default, the return of `DefaultAssay()` is 
#' used
#' @param rate_threshold Parameter to judge which the feature is expressed or
#' not in the cluster.
#' @param cluster_threshold Parameter to control expression specificity of
#' strong feature.
#'
#' @importFrom stats cor
#'
#' @returns Result data frame with DFC class.
#'
dfc_classify.Seurat <- function(
    data, dfc_res, assay = NULL,
    rate_threshold = 0.25, cluster_threshold = NULL
    ) {
  if(inherits(dfc_res,"dfc_models")) {
    stop("Please, input a DFC result object.")
  }
  dfc_res <- dfc_res$weights
  
  assay <- assay %||% DefaultAssay(data)
  smat <- GetAssayData(GetAssay(data,assay),'scale.data')
  if(length(smat)==0) {
    stop("Please, run Seurat::ScaleData")
  }
  cluster_label <- data$seurat_clusters
  useg <- intersect(dfc_res$feature,rownames(smat))
  smat <- t(smat[useg,])
  splited_smat <- split(as.data.frame(smat), f=cluster_label)
  posiRate <- sapply(splited_smat, function(x) {
    apply(x, 2, function(y) sum(y>min(y))/length(y))
    })
  posiCluster <- apply(posiRate>rate_threshold, 1, sum)
  if(is.null(cluster_threshold)) cluster_threshold <- floor(0.3*ncol(posiRate))
  dfc_class <- data.frame(
    weight = dfc_res$weight,
    posiCluster = posiCluster
    )
  dfc_class <- transform(
    dfc_class,
    class = ifelse(posiCluster == 0, "niche", posiCluster)
    )
  dfc_class <- transform(
    dfc_class,
    class = ifelse(posiCluster<cluster_threshold,
                   "strong", "weak")
    )
  return(dfc_class)
}

#' Specificity based DFC classification
#'
#' @param dfc_res A result of `dfc()`.
#' @param data The same one as when `dfc()` is run.
#' @param cluster_label A vector denoting the clusters of each cell
#' corresponding to the columns of data.
#' @param rate_threshold Parameter to judge which the feature is expressed or
#' not in the cluster.
#' @param cluster_threshold Parameter to control expression specificity of
#' strong feature.
#'
#' @importFrom stats cor
#'
#' @returns Result data frame with DFC class.
#'
dfc_classify.matrix <- function(
    data, dfc_res, cluster_label,
    rate_threshold = 0.25, cluster_threshold = NULL
    ) {
  if(inherits(dfc_res,"dfc_models")) {
    stop("Please, input a DFC result object.")
  }
  dfc_res <- dfc_res$weights
  
  useg <- intersect(dfc_res$feature,rownames(data))
  data <- t(data[useg,])
  splited_data <- split(as.data.frame(data), f=cluster_label)
  posiRate <- sapply(splited_data, function(x) {
    apply(x, 2, function(y) sum(y>min(y))/length(y))
  })
  posiCluster <- apply(posiRate>rate_threshold, 1, sum)
  if(is.null(cluster_threshold)) cluster_threshold <- floor(0.3*ncol(posiRate))
  dfc_class <- data.frame(
    weight = dfc_res$weight,
    posiCluster = posiCluster
  )
  dfc_class <- transform(
    dfc_class,
    class = ifelse(posiCluster == 0,
                   "niche", posiCluster)
  )
  dfc_class <- transform(
    dfc_class,
    class = ifelse(posiCluster < cluster_threshold,
                   "strong", "weak")
  )
  return(dfc_class)
}
