

#' Extract DFC
#'
#' @param data SeuratObject or gene x cell Matrix
#' @param ... ...
#'
#' @returns \item{Ridge}{The model of ridge regression to calculate penalry.}
#' \item{AdaLasso}{The model of Adaptive Lasso for DFC extraction.}
#' \item{weights}{A data frame of extracted features and the weights. When
#' return_Model set FALSE, only weigths are returned.}
#'
#' @export
#'
dfc <- function(data, ...) {
  UseMethod(generic = "dfc",object = data)
}

#' Specificity based DFC classification
#'
#' @param data The same one as when `dfc()` is run.'
#' @param ... ...
#'
#' @returns Result data frame with DFC class.
#'
#' @export
#'
dfc_classify <- function(data,...) {
  UseMethod(generic = "dfc_classify",object = data)
}

