

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
  if(!inherits(dfc_res,"dfc_models")) {
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
  if(!inherits(dfc_res,"dfc_models")) {
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
