#' legion package
#'
#' Package contains functions for multivariate time series forecasting
#'
#' \tabular{ll}{ Package: \tab smooth\cr Type: \tab Package\cr Date: \tab
#' 2016-01-27 - Inf\cr License: \tab GPL-2 \cr } The following functions are
#' included in the package:
#' \itemize{
#' \item \link[smooth]{ves} - Vector Exponential Smoothing.
#' \item \link[smooth]{ves} - Vector ETS-PIC model.
#' \item \link[smooth]{viss} - Multivariate occurrence ETS model.
#' }
#'
#' @name smooth
#' @docType package
#' @author Ivan Svetunkov
#' @author Kandrika Pritularga
#'
#' Maintainer: Ivan Svetunkov <ivan@svetunkov.ru>
#' @seealso \code{\link[greybox]{forecast}, \link[smooth]{es}, \link[smooth]{adam}}
#'
#' @template vssGeneralRef
#' @template vssKeywords
#'
#' @examples
#'
#' \dontrun{y <- cbind(rnorm(100,10,3),rnorm(100,10,3))
#'
#' ves(y,h=20,holdout=TRUE)
#'
#' @import zoo Rcpp
#' @importFrom nloptr nloptr
#' @importFrom graphics abline layout legend lines par points polygon
#' @importFrom stats AIC BIC cov end frequency is.ts median coef cor start time ts var simulate lm residuals
#' @importFrom utils packageVersion
#' @useDynLib legion
NULL



