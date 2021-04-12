#' @param loss Type of Loss Function used in optimization. \code{loss} can
#' be:
#' \itemize{
#' \item \code{"likelihood"} - which implies the maximisation of likelihood of
#' multivariate normal distribution (or log Normal if the multiplicative model
#' is constructed);
#' \item \code{"diagonal"} - similar to \code{"likelihood"}, but assumes that
#' covariances between the error terms are zero.
#' \item \code{"trace"} - the trace of the covariance matrix of errorrs.
#' The sum of variances is minimised in this case.
#' }
#' @param bounds What type of bounds to use in the model estimation. The first
#' letter can be used instead of the whole word. \code{"admissible"} means that the
#' model stability is ensured, while \code{"usual"} means that the all the parameters
#' are restricted by the (0, 1) region.
#' @param occurrence Defines type of occurrence model used. Can be:
#' \itemize{
#' \item \code{none}, meaning that the data should be considered as non-intermittent;
#' \item \code{fixed}, taking into account constant Bernoulli distribution of
#' demand occurrences;
#' \item \code{logistic}, based on logistic regression.
#' }
#' In this case, the ETS model inside the occurrence part will correspond to
#' \code{model} and \code{probability="dependent"}.
#' Alternatively, model estimated using \link[legion]{oves} function can be provided
#' here.
