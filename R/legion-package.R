#' Legion package
#'
#' Package contains functions for multivariate time series forecasting
#'
#' \tabular{ll}{ Package: \tab legion\cr Type: \tab Package\cr Date: \tab
#' 2021-02-18 - Inf\cr License: \tab GPL-2 \cr } The following functions are
#' included in the package:
#' \itemize{
#' \item \link[legion]{ves} - Vector Exponential Smoothing.
#' \item \link[legion]{vets} - Vector ETS-PIC model.
#' \item \link[legion]{oves} - Multivariate occurrence ETS model.
#' }
#'
#' @name legion
#' @docType package
#' @author Ivan Svetunkov
#' @author Kandrika Pritularga
#'
#' Maintainer: Ivan Svetunkov <ivan@svetunkov.ru>
#' @seealso \code{\link[generics]{forecast}, \link[smooth]{es}, \link[smooth]{adam}}
#'
#' @template vssGeneralRef
#' @template vssKeywords
#'
#' @examples
#'
#' \dontrun{y <- cbind(rnorm(100,10,3),rnorm(100,10,3))
#'
#' ves(y,h=20,holdout=TRUE)}
#'
#' @import zoo Rcpp
#' @importFrom nloptr nloptr
#' @importFrom graphics abline layout legend lines par points polygon
#' @importFrom stats AIC BIC cov median coef cor var lm residuals optim qt rnorm rt runif
#' @importFrom stats deltat end frequency is.ts start time ts
#' @importFrom utils packageVersion
#' @useDynLib legion
NULL



