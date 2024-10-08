% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/oves.R
\name{oves}
\alias{oves}
\title{Occurrence part of Vector State Space}
\usage{
oves(data, occurrence = c("logistic", "none", "fixed"), ic = c("AICc",
  "AIC", "BIC", "BICc"), h = 10, holdout = FALSE,
  probability = c("dependent", "independent"), model = "ANN",
  persistence = NULL, transition = NULL, phi = NULL, initial = NULL,
  initialSeason = NULL, xreg = NULL, ...)
}
\arguments{
\item{data}{The matrix with data, where series are in columns and
observations are in rows.}

\item{occurrence}{Type of method used in probability estimation. Can be
\code{"none"} - none, \code{"fixed"} - constant probability or
\code{"logistic"} - probability based on logit model.}

\item{ic}{Information criteria to use in case of model selection.}

\item{h}{Forecast horizon.}

\item{holdout}{If \code{TRUE}, holdout sample of size \code{h} is taken from
the end of the data.}

\item{probability}{Type of probability assumed in the model. If
\code{"dependent"}, then it is assumed that occurrence of one variable is
connected with the occurrence with another one. In case of \code{"independent"}
the occurrence of the variables is assumed to happen independent of each
other.}

\item{model}{Type of ETS model used for the estimation. Normally this should
be either \code{"ANN"} or \code{"MNN"}. If you assume that there are some
tendencies in occurrence, then you can use more complicated models. Model
selection is not yet available.}

\item{persistence}{Persistence matrix type. If \code{NULL}, then it is estimated.
See \link[legion]{ves} for the details.}

\item{transition}{Transition matrix type. If \code{NULL}, then it is estimated.
See \link[legion]{ves} for the details.}

\item{phi}{Damping parameter type. If \code{NULL}, then it is estimated.
See \link[legion]{ves} for the details.}

\item{initial}{Initial vector type. If \code{NULL}, then it is estimated.
See \link[legion]{ves} for the details.}

\item{initialSeason}{Type of the initial vector of seasonal components.
If \code{NULL}, then it is estimated. See \link[legion]{ves} for the details.}

\item{xreg}{Vector of matrix of exogenous variables, explaining some parts
of occurrence variable (probability).}

\item{...}{Other parameters. This is not needed for now.}
}
\value{
The object of class "oves" is returned. It contains following list of
values:

\itemize{
\item \code{model} - the type of the estimated ETS model;
\item \code{fitted} - fitted values of the constructed model;
\item \code{forecast} - forecast for \code{h} observations ahead;
\item \code{states} - values of states (currently level only);
\item \code{variance} - conditional variance of the forecast;
\item \code{logLik} - likelihood value for the model
\item \code{nParam} - number of parameters used in the model;
\item \code{residuals} - residuals of the model;
\item \code{data} - actual values of probabilities (zeros and ones).
\item \code{persistence} - the vector of smoothing parameters;
\item \code{initial} - initial values of the state vector;
\item \code{initialSeason} - the matrix of initials seasonal states;
\item \code{occurrence} - type of occurrence model used;
\item \code{probability} - type of probability used;
\item \code{issModel} - intermittent state-space model used for
calculations. Useful only in the case of \code{occurrence="l"} and
\code{probability="d"}.
}
}
\description{
Function calculates the probability for the occurrence part of vector state space model.
This is needed in order to forecast intermittent demand using other functions.
}
\details{
The function estimates probability of demand occurrence, using one of the VES
state-space models.
}
\examples{

Y <- cbind(c(rpois(25,0.1),rpois(25,0.5),rpois(25,1),rpois(25,5)),
           c(rpois(25,0.1),rpois(25,0.5),rpois(25,1),rpois(25,5)))

oves(Y, occurrence="l")
oves(Y, occurrence="l", probability="i")

}
\seealso{
\code{\link[smooth]{oes}, \link[legion]{ves}}
}
\author{
Ivan Svetunkov, \email{ivan@svetunkov.com}
}
\keyword{models}
\keyword{multivariate}
\keyword{nonlinear}
\keyword{regression}
\keyword{ts}
