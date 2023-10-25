

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
