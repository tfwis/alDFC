

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
#' @export
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
#' @export
#'
dfc.matrix <- function(
    data, target_clusters, cluster_label = NULL,
    gamma = 1, return_Model = FALSE, seed = NULL,
    lambda_penalty = "1se", lambda_weight = "1se",
    SIS = TRUE, min_feature = NULL, max_feature = NULL
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


