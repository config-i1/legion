% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/autovets.R, R/vets.R
\name{auto.vets}
\alias{auto.vets}
\alias{vets}
\title{Vector ETS-PIC model}
\usage{
auto.vets(data, model = "PPP", lags = c(frequency(data)),
  loss = c("likelihood", "diagonal", "trace"), ic = c("AICc", "AIC", "BIC",
  "BICc"), h = 10, holdout = FALSE, occurrence = c("none", "fixed",
  "logistic"), bounds = c("admissible", "usual", "none"), silent = TRUE,
  parallel = FALSE, ...)

vets(data, model = "PPP", lags = c(frequency(data)),
  parameters = c("level", "trend", "seasonal", "damped"),
  initials = c("seasonal"), components = c("none"),
  loss = c("likelihood", "diagonal", "trace"), ic = c("AICc", "AIC", "BIC",
  "BICc"), h = 10, holdout = FALSE, occurrence = c("none", "fixed",
  "logistic"), bounds = c("admissible", "usual", "none"), silent = TRUE,
  ...)
}
\arguments{
\item{data}{The matrix with the data, where series are in columns and
observations are in rows.}

\item{model}{The type of ETS model. Can consist of 3 or 4 chars: \code{ANN},
\code{AAN}, \code{AAdN}, \code{AAA}, \code{AAdA}, \code{MMdM} etc.
\code{PPP} means that the best pure model will be selected based on the chosen
information criteria type.
ATTENTION! ONLY PURE ADDITIVE AND PURE MULTIPLICATIVE MODELS ARE
AVAILABLE!
Pure multiplicative models are done as additive model applied to log(y).

Also \code{model} can accept a previously estimated VETS model and use all its
parameters.}

\item{lags}{The lags of the model. Needed for seasonal models.}

\item{loss}{Type of Loss Function used in optimization. \code{loss} can
be:
\itemize{
\item \code{"likelihood"} - which implies the maximisation of likelihood of
multivariate normal distribution (or log Normal if the multiplicative model
is constructed);
\item \code{"diagonal"} - similar to \code{"likelihood"}, but assumes that
covariances between the error terms are zero.
\item \code{"trace"} - the trace of the covariance matrix of errors.
The sum of variances is minimised in this case.
\item Provided by user as a custom function of \code{actual}, \code{fitted} and
\code{B}. Note that internally function transposes the data, so that \code{actual}
and \code{fitted} contain observations in columns and series in rows.
}
An example of the latter option is:
   \code{lossFunction <- function(actual,fitted,B){return(mean(abs(actual - fitted)))}}
followed by \code{loss=lossFunction}.}

\item{ic}{The information criterion used in the model selection procedure.}

\item{h}{Length of forecasting horizon.}

\item{holdout}{If \code{TRUE}, holdout sample of size \code{h} is taken from
the end of the data.}

\item{occurrence}{Defines type of occurrence model used. Can be:
\itemize{
\item \code{none}, meaning that the data should be considered as non-intermittent;
\item \code{fixed}, taking into account constant Bernoulli distribution of
demand occurrences;
\item \code{logistic}, based on logistic regression.
}
In this case, the ETS model inside the occurrence part will correspond to
\code{model} and \code{probability="dependent"}.
Alternatively, model estimated using \link[legion]{oves} function can be provided
here.}

\item{bounds}{What type of bounds to use in the model estimation. The first
letter can be used instead of the whole word. \code{"admissible"} means that the
model stability is ensured, while \code{"usual"} means that the all the parameters
are restricted by the (0, 1) region.}

\item{silent}{If \code{silent="none"}, then nothing is silent, everything is
printed out and drawn. \code{silent="all"} means that nothing is produced or
drawn (except for warnings). In case of \code{silent="graph"}, no graph is
produced. If \code{silent="legend"}, then legend of the graph is skipped.
And finally \code{silent="output"} means that nothing is printed out in the
console, but the graph is produced. \code{silent} also accepts \code{TRUE}
and \code{FALSE}. In this case \code{silent=TRUE} is equivalent to
\code{silent="all"}, while \code{silent=FALSE} is equivalent to
\code{silent="none"}. The parameter also accepts first letter of words ("n",
"a", "g", "l", "o").}

\item{parallel}{If TRUE, the estimation of ADAM models is done in parallel (used in \code{auto.adam} only).
If the number is provided (e.g. \code{parallel=41}), then the specified number of cores is set up.
WARNING! Packages \code{foreach} and either \code{doMC} (Linux and Mac only)
or \code{doParallel} are needed in order to run the function in parallel.}

\item{...}{Other non-documented parameters. For example \code{FI=TRUE} will
make the function also produce Fisher Information matrix, which then can be
used to calculated variances of smoothing parameters and initial states of
the model. The vector of initial parameter for the optimiser can be provided
here as the variable \code{B}. The upper bound for the optimiser is provided
via \code{ub}, while the lower one is \code{lb}. Also, the options for nloptr can be
passed here:
\itemize{
\item \code{maxeval=40*k} is the default number of iterations for both optimisers
used in the function (k is the number of parameters to estimate).
\item \code{algorithm1="NLOPT_LN_BOBYQA"} is the algorithm used in the first
optimiser, while \code{algorithm2="NLOPT_LN_NELDERMEAD"} is the second one.
\item \code{xtol_rel1=1e-8} is the relative tolerance in the first optimiser,
while \code{xtol_rel2=1e-6} is for the second one. All of this can be amended and
passed in ellipsis for finer tuning.
\item \code{print_level} - the level of output for the optimiser (0 by default).
If equal to 41, then the detailed results of the optimisation are returned.}}

\item{parameters}{The character vector, specifying, which of the parameters
should be common between time series. This includes smoothing parameters for
\code{"level"}, \code{"trend"}, \code{"seasonal"} components and \code{"damped"}
trend parameter. If \code{parameters="none"}, then all parameters are set to be
individual. An example is the model with all parameters being common:
\code{parameters=c("level","trend","seasonal","damped")}. The order is not important
and the first letters can be used instead of the full words as well.}

\item{initials}{The character vector, specifying, which of the initial values of
components should be common. This can be \code{"level"}, \code{"trend"} and / or
\code{"seasonal"}, setting initials of respective components to be common. This
can also be \code{"none"}, making the initials individual for all series. An
example is the model with only seasonal initials being common:
\code{initials="seasonal"}. The order is not important, and the first letters can
be used instead of the full words.}

\item{components}{The character vector, specifying, which of the components
components should be shared between time series. This can be \code{"level"},
\code{"trend"} and / or \code{"seasonal"}, setting respective components to be
shared. This can also be \code{"none"}, making them individual for all series.
The order is not important, and the first letters can be used instead of the full
words. Please, note that making components common automatically sets the
respective \code{initials} common as well.}
}
\value{
Object of class "legion" is returned. It contains the following list of
values:
\itemize{
\item \code{model} - The name of the fitted model;
\item \code{timeElapsed} - The time elapsed for the construction of the model;
\item \code{states} - The matrix of states with components in columns and time in rows;
\item \code{persistence} - The persistence matrix;
\item \code{transition} - The transition matrix;
\item \code{measurement} - The measurement matrix;
\item \code{phi} - The damping parameter value;
\item \code{B} - The vector of all the estimated coefficients;
\item \code{lagsAll} - The vector of the internal lags used in the model;
\item \code{nParam} - The number of estimated parameters;
\item \code{occurrence} - The occurrence model estimated with VETS;
\item \code{data} - The matrix with the original data;
\item \code{fitted} - The matrix of the fitted values;
\item \code{holdout} - The matrix with the holdout values (if \code{holdout=TRUE} in
the estimation);
\item \code{residuals} - The matrix of the residuals of the model;
\item \code{Sigma} - The covariance matrix of the errors (estimated with the correction
for the number of degrees of freedom);
\item \code{forecast} - The matrix of point forecasts;
\item \code{ICs} - The values of the information criteria;
\item \code{logLik} - The log-likelihood function;
\item \code{lossValue} - The value of the loss function;
\item \code{loss} - The type of the used loss function;
\item \code{lossFunction} - The loss function if the custom was used in the process;
\item \code{accuracy} - the values of the error measures. Currently not available.
\item \code{FI} - Fisher information if user asked for it using \code{FI=TRUE}.
}
}
\description{
Function constructs vector ETS model based on VETS-PIC taxonomy and returns
forecast, fitted values, errors and matrix of states along with other useful variables.
}
\details{
Function estimates vector ETS in the form of the Single Source of Error state space
model of the following type:

\deqn{
\mathbf{y}_{t} = \mathbf{o}_{t} (\mathbf{W} \mathbf{v}_{t-l} + \mathbf{x}_t
\mathbf{a}_{t-1} + \mathbf{\epsilon}_{t})
}{
y_{t} = o_{t} (W v_{t-l} + x_t a_{t-1} + \epsilon_{t})
}

\deqn{
\mathbf{v}_{t} = \mathbf{F} \mathbf{v}_{t-l} + \mathbf{G}
\mathbf{\epsilon}_{t}
}{
v_{t} = F v_{t-l} + G \epsilon_{t}
}

\deqn{\mathbf{a}_{t} = \mathbf{F_{X}} \mathbf{a}_{t-1} + \mathbf{G_{X}}
\mathbf{\epsilon}_{t} / \mathbf{x}_{t}}{a_{t} = F_{X} a_{t-1} + G_{X} \epsilon_{t}
/ x_{t}}

Where \eqn{y_{t}} is the vector of time series on observation \eqn{t}, \eqn{o_{t}}
is the vector of Bernoulli distributed random variable (in case of normal data it
becomes unit vector for all observations), \eqn{\mathbf{v}_{t}} is the matrix of
states and \eqn{l} is the matrix of lags, \eqn{\mathbf{x}_t} is the vector of
exogenous variables. \eqn{\mathbf{W}} is the measurement matrix, \eqn{\mathbf{F}}
is the transition matrix and \eqn{\mathbf{G}} is the persistence matrix.
Finally, \eqn{\epsilon_{t}} is the vector of error terms.

Conventionally we formulate values as:

\deqn{\mathbf{y}'_t = (y_{1,t}, y_{2,t}, \dots, y_{m,t})}{y_t = (y_{1,t}, y_{2,t},
\dots, y_{m,t}),}
where \eqn{m} is the number of series in the group.
\deqn{\mathbf{v}'_t = (v_{1,t}, v_{2,t}, \dots, v_{m,t})}{v'_t = (v_{1,t}, v_{2,t},
\dots, v_{m,t}),}
where \eqn{v_{i,t}} is vector of components for i-th time series.
\deqn{\mathbf{W}' = (w_{1}, \dots , 0;
\vdots , \ddots , \vdots;
0 , \vdots , w_{m})}{W' = (w_{1}, ... , 0;
... , ... , ...;
0 , ... , w_{m})} is matrix of measurement vectors.

The main idea of the function is in imposing restrictions on parameters / initials /
components of the model in order to capture the common dynamics between series.

In case of multiplicative model, instead of the vector y_t we use its logarithms.
As a result the multiplicative model is much easier to work with.

For some more information about the model and its implementation, see the
vignette: \code{vignette("vets","legion")}
}
\examples{

Y <- ts(cbind(rnorm(100,100,10),rnorm(100,75,8)),frequency=12)

# The simplest model applied to the data with the default values
vets(Y,model="ANN",h=10,holdout=TRUE)

# Multiplicative damped trend model with common parameters
# and initial seasonal indices
vets(Y,model="MMdM",h=10,holdout=TRUE,parameters=c("l","t","s","d"),
     initials="seasonal")

# Automatic selection of ETS components
\donttest{vets(Y, model="PPP", h=10, holdout=TRUE, initials="seasonal")}

}
\references{
\itemize{
\item de Silva A., Hyndman R.J. and Snyder, R.D. (2010). The vector
innovations structural time series framework: a simple approach to
multivariate forecasting. Statistical Modelling, 10 (4), pp.353-374
\item Hyndman, R.J., Koehler, A.B., Ord, J.K., and Snyder, R.D. (2008)
Forecasting with exponential smoothing: the state space approach,
Springer-Verlag.
\item Lütkepohl, H. (2005). New Introduction to Multiple Time Series
Analysis. New introduction to Multiple Time Series Analysis. Berlin,
Heidelberg: Springer Berlin Heidelberg.
\doi{10.1007/978-3-540-27752-1}
\item Svetunkov, I., Chen, H., & Boylan, J. E. (2023).
A new taxonomy for vector exponential smoothing and its application
to seasonal time series. European Journal of Operational Research,
304(3), 964–980. \doi{10.1016/j.ejor.2022.04.040}
}
}
\seealso{
\code{\link[legion]{ves}, \link[smooth]{es}, \link[smooth]{adam}}
}
\author{
Ivan Svetunkov, \email{ivan@svetunkov.com}
}
\keyword{models}
\keyword{multivariate}
\keyword{nonlinear}
\keyword{regression}
\keyword{ts}
