#' Specificity based DFC classification
#'
#' @param dfc_res A result of `dfc()`.
#' @param data The same one as when `dfc()` is run.
#' @param rate_threshold Parameter to judge which the feature is expressed or
#' not in the cluster.
#' @param cluster_threshold Parameter to control expression specificity of
#' strong feature.
#' @param ... ...
#'
#' @returns Result data frame with DFC class.
#'
#' @importFrom stats cor
#'
#' @export
#'
dfc_classify.matrix <- function(
    data, dfc_res, 
    rate_threshold = 0.25, cluster_threshold = NULL,
    eval_method = 'min', ...
) {
  if(!inherits(dfc_res,"dfc_models")) {
    stop("Please, input a DFC result object.")
  }
  dfc_weights <- dfc_res$weights
  cluster_label <- dfc_res$cluster_label
  target_label <- dfc_res$target_label
  useg <- intersect(dfc_weights$feature,rownames(data))
  data <- t(data[useg,])
  nCluster <- length(unique(cluster_label))
  if(is.null(cluster_threshold)) cluster_threshold <- floor(0.3*nCluster)
  eval_fun <- function(x) switch(eval_method,'min' = min(x),'zero' = 0)
  
  split_data <- split(as.data.frame(data), f=cluster_label)
  posiRate <- sapply(split_data, function(x) {
    apply(x, 2, function(y) sum(y>eval_fun(y))/length(y))
    })
  posiCluster <- apply(posiRate>rate_threshold, 1, sum)
  target_pn <- apply(data[target_label,], 2, function(y) {
    (sum(y>eval_fun(y))/length(y)) > rate_threshold
    })
  dfc_class <- data.frame(
    weight = dfc_weights$weight,
    posiCluster = posiCluster,
    expr_inTarget = target_pn
  )
  dfc_class <- transform(
    dfc_class,
    class = ifelse(posiCluster == 0, "niche",
                   ifelse(!target_pn,"weak",
                          ifelse(posiCluster < cluster_threshold,
                                 "strong", "weak")))
    )
  return(dfc_class)
}


#' Specificity based DFC classification
#'
#' @param data The same one as when `dfc()` is run.
#' @param dfc_res A result of `dfc()`.
#' @param assay The assay to use. In default, the return of `DefaultAssay()` is 
#' used
#' @param use.slot Please specify data slot to evaluate DFC genes.
#' @param ... ...
#'
#' @returns Result data frame with DFC class.
#' 
#' @importFrom stats cor
#' 
#' @export
#' 
dfc_classify.Seurat <- function(
    data, dfc_res, use.slot = 'scale.data', assay = NULL, ...
) {
  if(!inherits(dfc_res,"dfc_models")) {
    stop("Please, input a DFC result object.")
  }
  dfc_res_weight <- dfc_res$weights
  
  assay <- assay %||% DefaultAssay(data)
  smat <- GetAssayData(GetAssay(data,assay),slot = use.slot)
  if(length(smat)==0) {
    stop("Please, run Seurat::ScaleData")
  }
  if(!inherits(smat,"matrix")) {
    smat <- as.matrix(smat)
  }
  useg <- intersect(dfc_res_weight$feature,rownames(smat))
  smat <- as.matrix(smat[useg,])
  dfc_class <- dfc_classify(data = smat,dfc_res = dfc_res,...)
  return(dfc_class)
}


#' Specificity based DFC classification
#'
#' @param dfc_res A result of `dfc()`.
#' @param data The same one as when `dfc()` is run.
#' @param ... ...
#'
#' @returns Result data frame with DFC class.
#'
#' @importFrom stats cor
#'
#' @export
#'
dfc_classify.dgCMatrix <- function(
    data, dfc_res, ...
) {
  if(!inherits(dfc_res,"dfc_models")) {
    stop("Please, input a DFC result object.")
  }
  dfc_res_weights <- dfc_res$weights
  useg <- intersect(dfc_res$feature,rownames(data))
  data <- as.matrix(data[useg,])
  dfc_class <- dfc_classify(data = data,dfc_res = dfc_res, ...)
  return(dfc_class)
}

