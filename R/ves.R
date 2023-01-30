utils::globalVariables(c("nParamMax","nComponentsAll","nComponentsNonSeasonal","nSeries","modelIsSeasonal",
                         "obsInSample","obsAll","lagsModel","persistenceEstimate","persistenceType",
                         "persistenceValue","damped","dampedEstimate","dampedType","transitionType",
                         "initialEstimate","initialSeasonEstimate","initialSeasonValue","initialSeasonType",
                         "modelIsMultiplicative","matG","matW","B","ub","lb", "maxeval", "algorithm1",
                         "algorithm2", "xtol_rel1", "xtol_rel2", "Sigma","yFitted","PI","yDeltat",
                         "yFrequency","yStart","otObs","dataNames","seasonalType",
                         "CF","Etype","FI","ICs","Stype","Ttype","errors","h","holdout",
                         "initial","initialType","is.vsmooth.sim","lagsModel","lagsModelMax",
                         "level","matF","matVt","measures","nParam","normalizer","obsStates","ot",
                         "transition","transitionEstimate","yInSample","lossFunction",
                         "allowMultiplicative","modelDo","ICsAll",
                         "yClasses","yForecastStart","yInSampleIndex","yForecastIndex"));

#' Vector Exponential Smoothing in SSOE state space model
#'
#' Function constructs vector ETS model and returns forecast, fitted values, errors
#' and matrix of states along with other useful variables.
#'
#' Function estimates vector ETS in a form of the Single Source of Error state space
#' model of the following type:
#'
#' \deqn{
#' \mathbf{y}_{t} = (\mathbf{W} \mathbf{v}_{t-l} + \mathbf{x}_t
#' \mathbf{a}_{t-1} + \mathbf{\epsilon}_{t})
#' }{
#' y_{t} = (W v_{t-l} + x_t a_{t-1} + \epsilon_{t})
#' }
#'
#' \deqn{
#' \mathbf{v}_{t} = \mathbf{F} \mathbf{v}_{t-l} + \mathbf{G}
#' \mathbf{\epsilon}_{t}
#' }{
#' v_{t} = F v_{t-l} + G \epsilon_{t}
#' }
#'
#' \deqn{\mathbf{a}_{t} = \mathbf{F_{X}} \mathbf{a}_{t-1} + \mathbf{G_{X}}
#' \mathbf{\epsilon}_{t} / \mathbf{x}_{t}}{a_{t} = F_{X} a_{t-1} + G_{X} \epsilon_{t}
#' / x_{t}}
#'
#' Where \eqn{y_{t}} is the vector of time series on observation \eqn{t},
#' \eqn{\mathbf{v}_{t}} is the matrix of
#' states and \eqn{l} is the matrix of lags, \eqn{\mathbf{x}_t} is the vector of
#' exogenous variables. \eqn{\mathbf{W}} is the measurement matrix, \eqn{\mathbf{F}}
#' is the transition matrix and \eqn{\mathbf{G}} is the persistence matrix.
#' Finally, \eqn{\epsilon_{t}} is the vector of error terms.
#'
#' Conventionally we formulate values as:
#'
#' \deqn{\mathbf{y}'_t = (y_{1,t}, y_{2,t}, \dots, y_{m,t})}{y_t = (y_{1,t}, y_{2,t},
#' \dots, y_{m,t}),}
#' where \eqn{m} is the number of series in the group.
#' \deqn{\mathbf{v}'_t = (v_{1,t}, v_{2,t}, \dots, v_{m,t})}{v'_t = (v_{1,t}, v_{2,t},
#' \dots, v_{m,t}),}
#' where \eqn{v_{i,t}} is vector of components for i-th time series.
#' \deqn{\mathbf{W}' = (w_{1}, \dots , 0;
#' \vdots , \ddots , \vdots;
#' 0 , \vdots , w_{m})}{W' = (w_{1}, ... , 0;
#' ... , ... , ...;
#' 0 , ... , w_{m})} is matrix of measurement vectors.
#'
#' For the details on the additive model see Hyndman et al. (2008),
#' chapter 17.
#'
#' In case of multiplicative model, instead of the vector y_t we use its logarithms.
#' As a result the multiplicative model is much easier to work with.
#'
#' For some more information about the model and its implementation, see the
#' vignette: \code{vignette("ves","legion")}
#'
#' @template vssBasicParam
#' @template vssAdvancedParam
#' @template ssAuthor
#' @template vssKeywords
#'
#' @template vssGeneralRef
#'
#' @param model The type of ETS model. Can consist of 3 or 4 chars: \code{ANN},
#' \code{AAN}, \code{AAdN}, \code{AAA}, \code{AAdA}, \code{MMdM} etc.
#' \code{PPP} means that the best pure model will be selected based on the chosen
#' information criteria type.
#' ATTENTION! ONLY PURE ADDITIVE AND PURE MULTIPLICATIVE MODELS ARE
#' AVAILABLE!
#' Pure multiplicative models are done as additive model applied to log(y).
#'
#' Also \code{model} can accept a previously estimated VES model and use all its
#' parameters.
#' @param lags The lags of the model. Needed for seasonal models.
#' @param phi In cases of damped trend this parameter defines whether the \eqn{phi}
#' should be estimated separately for each series (\code{"individual"}) or for the whole
#' set (\code{"common"}). If vector or a value is provided here, then it is used by the
#' model.
#' @param initial Can be either character or a vector / matrix of initial states.
#' If it is character, then it can be \code{"individual"}, individual values of
#' the initial non-seasonal components are used, or \code{"common"}, meaning that
#' the initials for all the time series are set to be equal to the same value.
#' If vector of states is provided, then it is automatically transformed into
#' a matrix, assuming that these values are provided for the whole group.
#' @param initialSeason Can be either character or a vector / matrix of initial
#' states. Treated the same way as \code{initial}. This means that different time
#' series may share the same initial seasonal component.
# @param seasonal The type of seasonal component across the series. Can be
# \code{"individual"}, so that each series has its own component or \code{"common"},
# so that the component is shared across the series.
# @param weights The weights for the errors between the series with the common
# seasonal component. Ignored if \code{seasonal="individual"}.
#' @param persistence Persistence matrix \eqn{G}, containing smoothing
#' parameters. Can be:
#' \itemize{
#' \item \code{"independent"} - each series has its own smoothing parameters
#' and no interactions are modelled (all the other values in the matrix are set
#' to zero);
#' \item \code{"dependent"} - each series has its own smoothing parameters, but
#' interactions between the series are modelled (the whole matrix is estimated);
#' \item \code{"group"} each series has the same smoothing parameters for respective
#' components (the values of smoothing parameters are repeated, all the other values
#' in the matrix are set to zero).
#' \item \code{"seasonal"} - each component has its own smoothing parameter, except
#' for the seasonal one, which is common across the time series.
#' \item provided by user as a vector or as a matrix. The value is used by the model.
#' }
#' You can also use the first letter instead of writing the full word.
#' @param transition Transition matrix \eqn{F}. Can be:
#' \itemize{
#' \item \code{"independent"} - each series has its own preset transition matrix
#' and no interactions are modelled (all the other values in the matrix are set
#' to zero);
#' \item \code{"dependent"} - each series has its own transition matrix, but
#' interactions between the series are modelled (the whole matrix is estimated). The
#' estimated model behaves similar to VAR in this case;
#' \item \code{"group"} each series has the same transition matrix for respective
#' components (the values are repeated, all the other values in the matrix are set to
#' zero).
#' \item provided by user as a vector or as a matrix. The value is used by the model.
#' }
#' You can also use the first letter instead of writing the full word.
#' @param ...  Other non-documented parameters. For example \code{FI=TRUE} will
#' make the function also produce Fisher Information matrix, which then can be
#' used to calculated variances of smoothing parameters and initial states of
#' the model. The vector of initial parameter for the optimiser can be provided
#' here as the variable \code{B}. The upper bound for the optimiser is provided
#' via \code{ub}, while the lower one is \code{lb}. Also, the options for nloptr can be
#' passed here:
#' \itemize{
#' \item \code{maxeval=40*k} is the default number of iterations for both optimisers
#' used in the function (k is the number of parameters to estimate).
#' \item \code{algorithm1="NLOPT_LN_BOBYQA"} is the algorithm used in the first
#' optimiser, while \code{algorithm2="NLOPT_LN_NELDERMEAD"} is the second one.
#' \item \code{xtol_rel1=1e-8} is the relative tolerance in the first optimiser,
#' while \code{xtol_rel2=1e-6} is for the second one. All of this can be amended and
#' passed in ellipsis for finer tuning.
#' \item \code{print_level} - the level of output for the optimiser (0 by default).
#' If equal to 41, then the detailed results of the optimisation are returned.}
#' @return Object of class "legion" is returned. It contains the following list of
#' values:
#' \itemize{
#' \item \code{model} - The name of the fitted model;
#' \item \code{timeElapsed} - The time elapsed for the construction of the model;
#' \item \code{states} - The matrix of states with components in columns and time in rows;
#' \item \code{persistence} - The persistence matrix;
#' \item \code{transition} - The transition matrix;
#' \item \code{measurement} - The measurement matrix;
#' \item \code{phi} - The damping parameter value;
#' \item \code{lagsAll} - The vector of the internal lags used in the model;
#' \item \code{B} - The vector of all the estimated coefficients;
#' \item \code{initial} - The initial values of the non-seasonal components;
#' \item \code{initialSeason} - The initial values of the seasonal components;
#' \item \code{nParam} - The number of estimated parameters;
#' \item \code{occurrence} - The occurrence part of the model estimated with VES;
#' \item \code{data} - The matrix with the original data;
#' \item \code{fitted} - The matrix of the fitted values;
#' \item \code{holdout} - The matrix with the holdout values (if \code{holdout=TRUE} in
#' the estimation);
#' \item \code{residuals} - The matrix of the residuals of the model;
#' \item \code{Sigma} - The covariance matrix of the errors (estimated with the correction
#' for the number of degrees of freedom);
#' \item \code{forecast} - The matrix of point forecasts;
#' \item \code{ICs} - The values of the information criteria;
#' \item \code{logLik} - The log-likelihood function;
#' \item \code{lossValue} - The value of the loss function;
#' \item \code{loss} - The type of the used loss function;
#' \item \code{lossFunction} - The loss function if the custom was used in the process;
#' \item \code{accuracy} - the values of the error measures. Currently not available.
#' \item \code{FI} - Fisher information if user asked for it using \code{FI=TRUE}.
#' }
#' @seealso \code{\link[legion]{vets}, \link[smooth]{es}, \link[smooth]{adam}}
#'
#' @examples
#'
#' Y <- ts(cbind(rnorm(100,100,10),rnorm(100,75,8)),frequency=12)
#'
#' # The simplest model applied to the data with the default values
#' ves(Y,model="ANN",h=10,holdout=TRUE)
#'
#' # Damped trend model with the dependent persistence
#' ves(Y,model="AAdN",persistence="d",h=10,holdout=TRUE)
#'
#' # Multiplicative damped trend model with individual phi
#' \donttest{ves(Y,model="MMdM",persistence="i",h=10,holdout=TRUE,initialSeason="c")}
#'
#' # Automatic selection between pure models
#' \donttest{ves(Y,model="PPP",persistence="i",h=10,holdout=TRUE,initialSeason="c")}
#'
#' # Intermittent demand vector model
#' Y <- cbind(c(rpois(25,0.1),rpois(25,0.5),rpois(25,1),rpois(25,5)),
#'            c(rpois(25,0.1),rpois(25,0.5),rpois(25,1),rpois(25,5)))
#'
#' \donttest{ves(Y,model="MNN",h=10,holdout=TRUE,occurrence="l")}
#'
#' @importFrom smooth sowhat
#' @importFrom greybox measures
#' @export
ves <- function(data, model="PPP", lags=c(frequency(data)),
                persistence=c("common","individual","dependent"),
                transition=c("common","individual","dependent"), phi=c("common","individual"),
                initial=c("individual","common"), initialSeason=c("common","individual"),
                # seasonal=c("individual","common"), weights=rep(1/ncol(data),ncol(data)),
                loss=c("likelihood","diagonal","trace"),
                ic=c("AICc","AIC","BIC","BICc"), h=10, holdout=FALSE,
                occurrence=c("none","fixed","logistic"),
                bounds=c("admissible","usual","none"),
                silent=TRUE, ...){
    # Copyright (C) 2017 - Inf  Ivan Svetunkov

    # Start measuring the time of calculations
    startTime <- Sys.time();

    # If a previous model provided as a model, write down the variables
    if(any(is.legion(model))){
        persistence <- model$persistence;
        transition <- model$transition;
        phi <- model$phi;
        measurement <- model$measurement;
        initial <- model$initial;
        initialSeason <- model$initialSeason;
        # nParamOriginal <- model$nParam;
        # if(is.null(xreg)){
        #     xreg <- model$xreg;
        # }
        # initialX <- model$initialX;
        # persistenceX <- model$persistenceX;
        # transitionX <- model$transitionX;
        # if(any(c(persistenceX)!=0) | any((transitionX!=0)&(transitionX!=1))){
        #     updateX <- TRUE;
        # }
        model <- modelType(model);
    }
    # else{
    # nParamOriginal <- NULL;
    # }

    # Add all the variables in ellipsis to current environment
    list2env(list(...),environment());

    # Set up seasonal components as individual. Use vets() otherwise.
    seasonal <- "individual";
    weights <- rep(1/ncol(data),ncol(data));

    ##### Set environment for vssInput and make all the checks #####
    environment(vssInput) <- environment();
    vssInput("ves",ParentEnvironment=environment(),...);

    ##### Calculation of scale #####
    scalerVES <- function(distribution="dnorm", Etype, obsInSample, other=NULL,
                           errors, yFitted=NULL, normalizer=1, loss="likelihood"){
        if(loss=="likelihood"){
            scaleValue <- (errors / normalizer) %*% t(errors / normalizer) / obsInSample;
            return(scaleValue*normalizer^2);
        }
        else{
            scaleValue <- diag(rowSums(errors^2) / obsInSample);
            return(scaleValue);
        }
    }

    ##### Cost Function for VES #####
    CF <- function(B, loss="likelihood", Etype="A", Ttype="N", damped=FALSE,
                   nComponentsNonSeasonal, nComponentsAll, lagsModelMax){
        elements <- fillerVES(matVt, matF, matG, matW, B, Ttype, damped,
                              nComponentsNonSeasonal, nComponentsAll, lagsModelMax);

        # Check the bounds
        # Edit KFP: change symmetric to FALSE
        if(bounds=="a"){
            eigenValues <- eigen(elements$matF - elements$matG %*% elements$matW, only.values=TRUE, symmetric=FALSE)$values;
            if(max(abs(eigenValues)>(1 + 1E-50))){
                return(max(abs(eigenValues))*1E+100);
            }
        }

        # Fit the model
        fitting <- vFitterWrap(switch(Etype, "M"=log(yInSample), yInSample),
                               elements$matVt, elements$matF, elements$matW, elements$matG,
                               lagsModel, Etype, Ttype, Stype, ot);

        # Calculate the loss
        if(loss=="likelihood"){
            scaleValue <- scalerVES("dnorm", Etype, otObs, NULL,
                                     fitting$errors, NULL, normalizer=normalizer,
                                     loss="likelihood");

            cfRes <- -sum(switch(Etype,
                                 "A"=dmvnormInternal(fitting$errors, 0, scaleValue, log=TRUE),
                                 "M"=dmvnormInternal(fitting$errors, 0, scaleValue, log=TRUE)-
                                     colSums(log(yInSample))));
        }
        else if(loss=="GV"){
            scaleValue <- scalerVES("dnorm", Etype, otObs, NULL,
                                     fitting$errors, NULL, normalizer=normalizer,
                                     loss="likelihood");
            cfRes <- suppressWarnings(log(det(scaleValue)) + nSeries * log(normalizer^2));
        }
        else if(loss=="diagonal"){
            scaleValue <- scalerVES("dnorm", Etype, otObs, NULL,
                                     fitting$errors, NULL, normalizer=normalizer,
                                     loss="diagonal");
            cfRes <- -sum(switch(Etype,
                                 "A"=dmvnormInternal(fitting$errors, 0, scaleValue, log=TRUE),
                                 "M"=dmvnormInternal(fitting$errors, 0, scaleValue, log=TRUE)-
                                     colSums(log(yInSample))));
        }
        # Loss for the oves model
        else if(loss=="occurrence"){
            cfRes <- -sum(log(fitting$yfit[yInSample==1]));
        }
        # Custom loss
        else if(loss=="custom"){
            cfRes <- switch(Etype,
                            "A"=lossFunction(actual=yInSample,fitted=fitting$yfit,B=B),
                            "M"=lossFunction(actual=yInSample,fitted=exp(fitting$yfit),B=B));
        }
        else{
            cfRes <- sum(rowSums(fitting$errors^2)) / obsInSample;
        }

        if(is.nan(cfRes) | is.na(cfRes) | is.infinite(cfRes)){
            cfRes <- 1e+100;
        }

        return(cfRes);
    }

    ##### LogLik for VES #####
    logLikVES <- function(B, loss="likelihood", Etype="A"){
        if(any(loss==c("likelihood","diagonal"))){
            return(-CF(B, loss=loss, Etype=Etype, Ttype=Ttype, damped=damped,
                       nComponentsNonSeasonal, nComponentsAll, lagsModelMax));
        }
        else{
            return(-CF(B, loss="likelihood", Etype=Etype, Ttype=Ttype, damped=damped,
                       nComponentsNonSeasonal, nComponentsAll, lagsModelMax));
        }
    }

    ##### IC values for VETS #####
    ICsVES <- function(B, logLikVESValue, nSeries, nParamPerSeries, obsInSample){
        ICs <- setNames(vector("numeric",4),c("AIC","AICc","BIC","BICc"));
        nParamAll <- nparam(logLikVESValue);

        # AIC
        ICs[1] <- AIC(logLikVESValue);
        # AICc
        if(obsInSample - (nParamPerSeries + nSeries + 1) <=0){
            ICs[2] <- Inf;
        }
        else{
            ICs[2] <- -2*logLikVESValue + 2*(obsInSample*nParamAll /
                                                  (obsInSample - (nParamPerSeries + nSeries + 1)));
        }
        # BIC
        ICs[3] <- BIC(logLikVESValue);
        # BICc
        # All the estimated parameters (would differ depending on loss)
        if(obsInSample - (nParamPerSeries + nSeries + 1) <=0){
            ICs[4] <- Inf;
        }
        else{
            ICs[4] <- -2*logLikVESValue + log(obsInSample)*(obsInSample*nParamAll /
                                                                 (obsInSample - (nParamPerSeries + nSeries + 1)));
        }
        return(ICs);
    }

    ##### B values for estimation #####
    # Function constructs default bounds where B values should lie
    initialiserVES <- function(Ttype,Stype,lagsModelMax,nComponentsAll,nComponentsNonSeasonal,nSeries){
        B <- NA;
        BLower <- NA;
        BUpper <- NA;
        BNames <- NA;

        if(seasonalType=="i"){
            #### Individual seasonality ####
            ### Persistence matrix
            # Edit KFP: change the parameter space, make all persistence to have mean-reverting behaviour
            # by forcing the off-diagonals to be negative so that the MA coefficient would be positive
            # potentially works for VES(ANN) only for now - dependent persistence matrix
            if(persistenceEstimate){
                if(persistenceType=="c"){
                    persistenceLength <- nComponentsAll;
                    if(bounds=="u"){
                        BLower <- c(BLower,rep(0,persistenceLength));
                        BUpper <- c(BUpper,rep(1,persistenceLength));
                    }
                    else{
                        BLower <- c(BLower,rep(-5,persistenceLength));
                        BUpper <- c(BUpper,rep(5,persistenceLength));
                    }
                }
                else if(persistenceType=="i"){
                    persistenceLength <- nComponentsAll*nSeries;
                    if(bounds=="u"){
                        BLower <- c(BLower,rep(0,persistenceLength));
                        BUpper <- c(BUpper,rep(1,persistenceLength));
                    }
                    else{
                        BLower <- c(BLower,rep(-5,persistenceLength));
                        BUpper <- c(BUpper,rep(5,persistenceLength));
                    }
                }
                # Edited by KFP
                else if(persistenceType=="d"){
                    persistenceLength <- nComponentsAll*nSeries^2;
                    if(bounds=="u"){
                        BLower <- c(BLower,rep(0,persistenceLength));
                        BUpper <- c(BUpper,rep(1,persistenceLength));
                    }
                    else{
                        BLower <- c(BLower,rep(-5,persistenceLength));
                        BUpper <- c(BUpper,rep(5,persistenceLength));
                    }
                }
                else if(persistenceType=="s"){
                    persistenceLength <- (nComponentsAll-1)*nSeries+1;
                    if(bounds=="u"){
                        BLower <- c(BLower,rep(0,persistenceLength));
                        BUpper <- c(BUpper,rep(1,persistenceLength));
                    }
                    else{
                        BLower <- c(BLower,rep(-5,persistenceLength));
                        BUpper <- c(BUpper,rep(5,persistenceLength));
                    }
                }
                B <- c(B,rep(1/persistenceLength,persistenceLength));

                BNames <- c(BNames,paste0("Persistence",c(1:persistenceLength)));
            }

            ### Damping parameter
            if(dampedEstimate){
                if(dampedType=="c"){
                    dampedLength <- 1;
                }
                else if(dampedType=="i"){
                    dampedLength <- nSeries;
                }
                B <- c(B,rep(0.95,dampedLength));
                BLower <- c(BLower,rep(0,dampedLength));
                BUpper <- c(BUpper,rep(1,dampedLength));
                BNames <- c(BNames,paste0("phi",c(1:dampedLength)));
            }

            ### Transition matrix
            if(transitionEstimate){
                if(transitionType=="d"){
                    transitionLength <- ((nSeries-1)*nComponentsAll^2)*nSeries;
                }
                B <- c(B,rep(0.1,transitionLength));
                BLower <- c(BLower,rep(-1,transitionLength));
                BUpper <- c(BUpper,rep(1,transitionLength));
                BNames <- c(BNames,paste0("transition",c(1:transitionLength)));
            }

            ### Vector of initials
            if(initialEstimate){
                if(initialType=="c"){
                    initialLength <- nComponentsNonSeasonal;
                }
                else{
                    initialLength <- nComponentsNonSeasonal*nSeries;
                }
                B <- c(B,initialValue);
                BNames <- c(BNames,paste0("initial",c(1:initialLength)));
                BLower <- c(BLower,rep(-Inf,initialLength));
                BUpper <- c(BUpper,rep(Inf,initialLength));
            }

            ### Vector of initial seasonals
            if(modelIsSeasonal && initialSeasonEstimate){
                if(initialSeasonType=="c"){
                    initialSeasonLength <- lagsModelMax;
                }
                else{
                    initialSeasonLength <- lagsModelMax*nSeries;
                }
                B <- c(B,initialSeasonValue);
                BNames <- c(BNames,paste0("initialSeason",c(1:initialSeasonLength)));
                # if(Stype=="A"){
                BLower <- c(BLower,rep(-Inf,initialSeasonLength));
                BUpper <- c(BUpper,rep(Inf,initialSeasonLength));
                # }
                # else{
                #     BLower <- c(BLower,rep(-0.0001,initialSeasonLength));
                #     BUpper <- c(BUpper,rep(20,initialSeasonLength));
                # }
            }
        }
        else{
            #### Common seasonality ####
            ### Persistence matrix
            if(persistenceEstimate){
                if(persistenceType=="c"){
                    persistenceLength <- nComponentsAll;
                }
                else if(persistenceType=="i"){
                    persistenceLength <- nComponentsNonSeasonal*nSeries+nSeries;
                }
                else if(persistenceType=="d"){
                    persistenceLength <- nComponentsNonSeasonal*nSeries^2+nSeries;
                }
                else if(persistenceType=="s"){
                    persistenceLength <- nComponentsNonSeasonal*nSeries+1;
                }
                B <- c(B,rep(0.1,persistenceLength));
                if(bounds=="u"){
                    BLower <- c(BLower,rep(0,persistenceLength));
                    BUpper <- c(BUpper,rep(1,persistenceLength));
                }
                else{
                    BLower <- c(BLower,rep(-5,persistenceLength));
                    BUpper <- c(BUpper,rep(5,persistenceLength));
                }
                BNames <- c(BNames,paste0("Persistence",c(1:persistenceLength)));
            }

            ### Damping parameter
            if(dampedEstimate){
                if(dampedType=="c"){
                    dampedLength <- 1;
                }
                else if(dampedType=="i"){
                    dampedLength <- nSeries;
                }
                B <- c(B,rep(0.95,dampedLength));
                BLower <- c(BLower,rep(0,dampedLength));
                BUpper <- c(BUpper,rep(1,dampedLength));
                BNames <- c(BNames,paste0("phi",c(1:dampedLength)));
            }

            ### Transition matrix
            if(transitionEstimate){
                if(transitionType=="d"){
                    transitionLength <- ((nSeries-1)*nComponentsNonSeasonal^2)*nSeries;
                }
                B <- c(B,rep(0.1,transitionLength));
                BLower <- c(BLower,rep(-1,transitionLength));
                BUpper <- c(BUpper,rep(1,transitionLength));
                BNames <- c(BNames,paste0("transition",c(1:transitionLength)));
            }

            ### Vector of initials
            if(initialEstimate){
                if(initialType=="c"){
                    initialLength <- nComponentsNonSeasonal;
                }
                else{
                    initialLength <- nComponentsNonSeasonal*nSeries;
                }
                B <- c(B,initialValue);
                BNames <- c(BNames,paste0("initial",c(1:initialLength)));
                BLower <- c(BLower,rep(-Inf,initialLength));
                BUpper <- c(BUpper,rep(Inf,initialLength));
            }

            ### Vector of initial seasonals
            if(initialSeasonEstimate){
                initialSeasonLength <- lagsModelMax;
                B <- c(B,initialSeasonValue);
                BNames <- c(BNames,paste0("initialSeason",c(1:initialSeasonLength)));
                BLower <- c(BLower,rep(-Inf,initialSeasonLength));
                BUpper <- c(BUpper,rep(Inf,initialSeasonLength));
            }
        }

        B <- B[!is.na(B)];
        BLower <- BLower[!is.na(BLower)];
        BUpper <- BUpper[!is.na(BUpper)];
        BNames <- BNames[!is.na(BNames)];

        return(list(B=B,BLower=BLower,BUpper=BUpper,BNames=BNames));
    }

    ##### Basic VES architector #####
    ### This function will accept Etype, Ttype, Stype and damped and would return:
    # nComponentsNonSeasonal, nComponentsAll, lagsModelMax, modelIsSeasonal, obsStates
    # This is needed for model selection
    architectorVES <- function(Etype, Ttype, Stype, damped, nSeries, lags){
        # Binaries for trend and seasonal
        modelIsTrendy <- Ttype!="N";
        modelIsSeasonal <- Stype!="N";

        lagsModelMax <- max(lags) * modelIsSeasonal + 1 * (!modelIsSeasonal);

        # Define the number of rows that should be in the matVt
        obsStates <- max(obsAll + lagsModelMax, obsInSample + 2*lagsModelMax);

        nComponentsNonSeasonal <- 1 + (Ttype!="N")*1;
        nComponentsAll <- nComponentsNonSeasonal + modelIsSeasonal*1;

        return(list(modelIsTrendy=modelIsTrendy, modelIsSeasonal=modelIsSeasonal, lagsModelMax=lagsModelMax,
                    nComponentsNonSeasonal=nComponentsNonSeasonal, nComponentsAll=nComponentsAll));
    }

    ##### Basic matrices creator #####
    # This thing returns matVt, matF, matG, matW, dampedValue, initialValue
    # and initialSeasonValue if they are not provided + lagsModel
    creatorVES <- function(...){
        # ellipsis <- list(...);
        # ParentEnvironment <- ellipsis[['ParentEnvironment']];

        ### Persistence matrix
        matG <- switch(seasonalType,
                       "i" =  matrix(0,nSeries*nComponentsAll,nSeries),
                       "c" = matrix(0,nSeries*nComponentsNonSeasonal+1,nSeries));
        if(!persistenceEstimate){
            matG <- persistenceValue;
        }

        ### Damping parameter
        if(!damped){
            dampedValue <- matrix(1,nSeries,1);
        }

        ### Transition matrix
        if(any(transitionType==c("c","i","d"))){
            if(Ttype=="N"){
                transitionValue <- matrix(1,1,1);
            }
            else if(Ttype!="N"){
                transitionValue <- matrix(c(1,0,dampedValue[1],dampedValue[1]),2,2);
            }
            if(Stype!="N"){
                transitionValue <- cbind(transitionValue,rep(0,nComponentsNonSeasonal));
                transitionValue <- rbind(transitionValue,c(rep(0,nComponentsNonSeasonal),1));
            }
            transitionValue <- matrix(transitionValue,nComponentsAll,nComponentsAll);
            transitionBuffer <- diag(nSeries*nComponentsAll);
            for(i in 1:nSeries){
                transitionBuffer[c(1:nComponentsAll)+nComponentsAll*(i-1),
                                 c(1:nComponentsAll)+nComponentsAll*(i-1)] <- transitionValue;
            }
            if(any(transitionType==c("i","d")) & damped){
                for(i in 1:nSeries){
                    transitionBuffer[c(1:nComponentsNonSeasonal)+nComponentsAll*(i-1),
                                     nComponentsNonSeasonal+nComponentsAll*(i-1)] <- dampedValue[i];
                }
            }
            transitionValue <- transitionBuffer;
        }
        if(transitionType=="d"){
            # Fill in the other values of F with some values
            for(i in 1:nSeries){
                transitionValue[c(1:nComponentsAll)+nComponentsAll*(i-1),
                                setdiff(c(1:nSeries*nComponentsAll),c(1:nComponentsAll)+nComponentsAll*(i-1))] <- 0.1;
            }
        }
        matF <- switch(seasonalType,
                       "i"=transitionValue,
                       "c"=rbind(cbind(transitionValue[-(c(1:nSeries)*nComponentsAll),
                                                       -(c(1:nSeries)*nComponentsAll)],
                                       0),
                                 c(transitionValue[nComponentsAll*nSeries,
                                                   -(c(1:nSeries)*nComponentsAll)],1)));

        ### Measurement matrix
        if(seasonalType=="i"){
            matW <- matrix(0,nSeries,nSeries*nComponentsAll);
            for(i in 1:nSeries){
                matW[i,c(1:nComponentsAll)+nComponentsAll*(i-1)] <- 1;
            }
            if(damped){
                for(i in 1:nSeries){
                    matW[i,nComponentsNonSeasonal+nComponentsAll*(i-1)] <- dampedValue[i];
                }
            }
        }
        else{
            matW <- matrix(0,nSeries,nSeries*nComponentsNonSeasonal+1);
            for(i in 1:nSeries){
                matW[i,c(1:nComponentsNonSeasonal)+nComponentsNonSeasonal*(i-1)] <- 1;
            }
            matW[,nSeries*nComponentsNonSeasonal+1] <- 1;
            if(damped){
                for(i in 1:nSeries){
                    matW[i,nComponentsNonSeasonal+nComponentsNonSeasonal*(i-1)] <- dampedValue[i];
                }
            }
        }

        ### Vector of states
        statesNames <- "level";
        if(Ttype!="N"){
            statesNames <- c(statesNames,"trend");
        }
        if(Stype!="N"){
            statesNames <- c(statesNames,"seasonal");
        }
        matVt <- matrix(NA, nComponentsAll*nSeries, obsStates,
                        dimnames=list(paste0(rep(dataNames,each=nComponentsAll),
                                             "_",statesNames),NULL));
        if(seasonalType=="c"){
            matVt <- rbind(matVt[-(c(1:nSeries)*nComponentsAll),,drop=F],
                           matVt[nComponentsAll,,drop=F]);
            rownames(matVt)[nComponentsNonSeasonal*nSeries+1] <- "seasonal";
        }

        ## Deal with non-seasonal part of the vector of states
        if(!initialEstimate){
            if(seasonalType=="i"){
                initialPlaces <- nComponentsAll*(c(1:nSeries)-1)+1;
                if(Ttype!="N"){
                    initialPlaces <- c(initialPlaces,nComponentsAll*(c(1:nSeries)-1)+2);
                    initialPlaces <- sort(initialPlaces);
                }
            }
            else{
                initialPlaces <- c(1:(nSeries*nComponentsNonSeasonal));
            }
            matVt[initialPlaces,1:lagsModelMax] <- rep(initialValue,lagsModelMax);
        }
        else{
            XValues <- rbind(rep(1,obsInSample),c(1:obsInSample));
            initialValue <- switch(Etype, "M"=log(yInSample), yInSample) %*% t(XValues) %*% solve(XValues %*% t(XValues));
            if(Etype=="L"){
                initialValue[,1] <- (initialValue[,1] - 0.5) * 20;
            }

            #### !!! This is a temporary fix for log(0) on intermittent demand ####
            if(any(is.nan(initialValue))){
                initialValue[,1] <- rowMeans(yInSample);
                initialValue[,2] <- switch(Ttype, "M"=1, 0);
            }

            if(Ttype=="N"){
                initialValue <- matrix(initialValue[,-2],nSeries,1);
            }
            if(initialType=="c"){
                initialValue <- matrix(colMeans(initialValue),nComponentsNonSeasonal,1);
            }
            else{
                initialValue <- matrix(as.vector(t(initialValue)),nComponentsNonSeasonal * nSeries,1);
            }
        }

        ## Deal with seasonal part of the vector of states
        if(modelIsSeasonal){
            if(initialSeasonType=="p"){
                if(seasonalType=="i"){
                    initialPlaces <- nComponentsAll*(c(1:nSeries)-1)+nComponentsAll;
                }
                else{
                    initialPlaces <- nSeries*nComponentsNonSeasonal+1;
                }
                matVt[initialPlaces,1:lagsModelMax] <- initialSeasonValue;
            }
            else{
                # Matrix of dummies for seasons
                XValues <- matrix(rep(diag(lagsModelMax),ceiling(obsInSample/lagsModelMax)),lagsModelMax)[,1:obsInSample];
                # if(Stype=="A"){
                initialSeasonValue <- (switch(Etype, "M"=log(yInSample), yInSample)-
                                           rowMeans(switch(Etype, "M"=log(yInSample), yInSample))) %*%
                    t(XValues) %*% solve(XValues %*% t(XValues));
                # }
                # else{
                #     initialSeasonValue <- (yInSample-rowMeans(yInSample)) %*% t(XValues) %*% solve(XValues %*% t(XValues));
                # }
                if(initialSeasonType=="c" || seasonalType=="c"){
                    initialSeasonValue <- matrix(colMeans(initialSeasonValue),1,lagsModelMax);
                }
                else{
                    initialSeasonValue <- matrix(as.vector(t(initialSeasonValue)),nSeries,lagsModelMax);
                }
            }
        }

        ### lagsModel
        if(seasonalType=="i"){
            lagsModel <- rep(1,nComponentsAll);
            if(modelIsSeasonal){
                lagsModel[nComponentsAll] <- lagsModelMax;
            }
            lagsModel <- matrix(lagsModel,nSeries*nComponentsAll,1);
        }
        else{
            lagsModel <- matrix(c(rep(1,nSeries*nComponentsNonSeasonal),lagsModelMax),
                                nSeries*nComponentsNonSeasonal+1,1);
        }

        return(list(matVt=matVt, matF=matF, matG=matG, matW=matW, dampedValue=dampedValue,
                    initialValue=initialValue, initialSeasonValue=initialSeasonValue, lagsModel=lagsModel));
    }

    ##### Basic matrices filler #####
    # This thing fills in matVt, matF, matG and matW with values from B and returns the corrected values
    fillerVES <- function(matVt,matF,matG,matW,B,Ttype,damped,
                          nComponentsNonSeasonal,nComponentsAll,lagsModelMax){
        nCoefficients <- 0;
        ##### Individual seasonality #####
        if(seasonalType=="i"){
            ### Persistence matrix
            if(persistenceEstimate){
                persistenceBuffer <- matrix(0,nSeries*nComponentsAll,nSeries);
                # Grouped values
                if(persistenceType=="c"){
                    persistenceValue <- B[1:nComponentsAll];
                    nCoefficients <- nComponentsAll;
                    for(i in 1:nSeries){
                        persistenceBuffer[1:nComponentsAll+nComponentsAll*(i-1),i] <- persistenceValue;
                    }
                    persistenceValue <- persistenceBuffer;
                }
                # Independent values
                else if(persistenceType=="i"){
                    persistenceValue <- B[1:(nComponentsAll*nSeries)];
                    nCoefficients <- nComponentsAll*nSeries;
                    for(i in 1:nSeries){
                        persistenceBuffer[1:nComponentsAll+nComponentsAll*(i-1),
                                          i] <- persistenceValue[1:nComponentsAll+nComponentsAll*(i-1)];
                    }
                    persistenceValue <- persistenceBuffer;
                }
                # Dependent values
                else if(persistenceType=="d"){
                    persistenceValue <- B[1:(nComponentsAll*nSeries^2)];
                    nCoefficients <- nComponentsAll*nSeries^2;
                }
                # Grouped seasonal values
                else if(persistenceType=="s"){
                    persistenceValue <- B[1:((nComponentsAll-1)*nSeries+1)];
                    persistenceSeasonal <- persistenceValue[length(persistenceValue)];
                    nCoefficients <- ((nComponentsAll-1)*nSeries+1);
                    for(i in 1:nSeries){
                        persistenceBuffer[1:(nComponentsAll-1)+nComponentsAll*(i-1),
                                          i] <- persistenceValue[1:(nComponentsAll-1)+(nComponentsAll-1)*(i-1)];
                        persistenceBuffer[nComponentsAll+nComponentsAll*(i-1),i] <- persistenceSeasonal;
                    }
                    persistenceValue <- persistenceBuffer;
                }
                matG[,] <- persistenceValue;
            }

            ### Damping parameter
            if(modelIsTrendy && damped){
                if(dampedType=="c"){
                    dampedValue <- matrix(B[nCoefficients+1],nSeries,1);
                    nCoefficients <- nCoefficients + 1;
                }
                else if(dampedType=="i"){
                    dampedValue <- matrix(B[nCoefficients+(1:nSeries)],nSeries,1);
                    nCoefficients <- nCoefficients + nSeries;
                }
            }

            ### Transition matrix
            if(any(transitionType==c("i","d","c")) && modelIsTrendy && damped){
                for(i in 1:nSeries){
                    matF[c(1:nComponentsNonSeasonal)+nComponentsAll*(i-1),
                         nComponentsNonSeasonal+nComponentsAll*(i-1)] <- dampedValue[i];
                }
            }
            if(transitionType=="d"){
                # Fill in the other values of F with some values
                nCoefficientsBuffer <- (nSeries-1)*nComponentsAll^2;

                for(i in 1:nSeries){
                    matF[c(1:nComponentsAll)+nComponentsAll*(i-1),
                         setdiff(c(1:(nSeries*nComponentsAll)),
                                 c(1:nComponentsAll)+nComponentsAll*(i-1))] <- B[nCoefficients+c(1:nCoefficientsBuffer)];
                    nCoefficients <- nCoefficients + nCoefficientsBuffer;
                }
            }

            ### Measurement matrix
            # Needs to be filled in with dampedValue even if dampedValue has been provided by a user
            if(modelIsTrendy && damped){
                for(i in 1:nSeries){
                    matW[i,nComponentsNonSeasonal+nComponentsAll*(i-1)] <- dampedValue[i];
                }
            }

            ### Vector of states
            ## Deal with non-seasonal part of the vector of states
            if(initialEstimate){
                initialPlaces <- nComponentsAll*(c(1:nSeries)-1)+1;
                if(Ttype!="N"){
                    initialPlaces <- c(initialPlaces,nComponentsAll*(c(1:nSeries)-1)+2);
                    initialPlaces <- sort(initialPlaces);
                }
                if(initialType=="i"){
                    initialValue <- matrix(B[nCoefficients+c(1:(nComponentsNonSeasonal*nSeries))],
                                           nComponentsNonSeasonal * nSeries,1);
                    nCoefficients <- nCoefficients + nComponentsNonSeasonal*nSeries;
                }
                else if(initialType=="c"){
                    initialValue <- matrix(B[nCoefficients+c(1:nComponentsNonSeasonal)],
                                           nComponentsNonSeasonal * nSeries,1);
                    nCoefficients <- nCoefficients + nComponentsNonSeasonal;
                }
                matVt[initialPlaces,1:lagsModelMax] <- rep(initialValue,lagsModelMax);
            }

            ## Deal with seasonal part of the vector of states
            if(modelIsSeasonal & initialSeasonEstimate){
                initialPlaces <- nComponentsAll*(c(1:nSeries)-1)+nComponentsAll;
                if(initialSeasonType=="i"){
                    matVt[initialPlaces,1:lagsModelMax] <- matrix(B[nCoefficients+c(1:(nSeries*lagsModelMax))],
                                                                  nSeries,lagsModelMax,byrow=TRUE);
                    nCoefficients <- nCoefficients + nSeries*lagsModelMax;
                }
                else if(initialSeasonType=="c"){
                    matVt[initialPlaces,1:lagsModelMax] <- matrix(B[nCoefficients+c(1:lagsModelMax)],nSeries,lagsModelMax,byrow=TRUE);
                    nCoefficients <- nCoefficients + lagsModelMax;
                }
            }
        }
        ##### Common seasonality #####
        else{
            ### Persistence matrix
            if(persistenceEstimate){
                persistenceBuffer <- matrix(0,nSeries*nComponentsNonSeasonal+1,nSeries);
                # Grouped values
                if(persistenceType=="c"){
                    persistenceValue <- B[1:nComponentsAll];
                    nCoefficients <- nComponentsAll;
                    for(i in 1:nSeries){
                        persistenceBuffer[1:nComponentsNonSeasonal+nComponentsNonSeasonal*(i-1),
                                          i] <- persistenceValue[1:nComponentsNonSeasonal];
                    }
                    persistenceBuffer[nSeries*nComponentsNonSeasonal+1,] <- weights*persistenceValue[nComponentsAll];
                    persistenceValue <- persistenceBuffer;
                }
                # Independent values
                else if(persistenceType=="i"){
                    persistenceValue <- B[1:(nComponentsNonSeasonal*nSeries+nSeries)];
                    nCoefficients <- nComponentsNonSeasonal*nSeries+nSeries;
                    for(i in 1:nSeries){
                        persistenceBuffer[1:nComponentsNonSeasonal+nComponentsNonSeasonal*(i-1),
                                          i] <- persistenceValue[1:nComponentsNonSeasonal+nComponentsNonSeasonal*(i-1)];
                    }
                    persistenceBuffer[nSeries*nComponentsNonSeasonal+1,
                                      ] <- weights*persistenceValue[nComponentsNonSeasonal*nSeries+c(1:nSeries)];
                    persistenceValue <- persistenceBuffer;
                }
                # Dependent values
                else if(persistenceType=="d"){
                    persistenceValue <- B[1:(nSeries^2*nComponentsNonSeasonal+nSeries)];
                    nCoefficients <- nSeries^2*nComponentsNonSeasonal+nSeries;
                }
                # Grouped seasonal values
                else if(persistenceType=="s"){
                    persistenceValue <- B[1:(nComponentsNonSeasonal*nSeries+1)];
                    nCoefficients <- nComponentsNonSeasonal*nSeries+1;
                    for(i in 1:nSeries){
                        persistenceBuffer[1:nComponentsNonSeasonal+nComponentsNonSeasonal*(i-1),
                                          i] <- persistenceValue[1:nComponentsNonSeasonal+nComponentsNonSeasonal*(i-1)];
                    }
                    persistenceBuffer[nSeries*nComponentsNonSeasonal+1,
                                      ] <- weights*persistenceValue[nComponentsNonSeasonal*nSeries+1];
                    persistenceValue <- persistenceBuffer;
                }
                matG[,] <- persistenceValue;
            }

            ### Damping parameter
            if(damped){
                if(dampedType=="c"){
                    dampedValue <- matrix(B[nCoefficients+1],nSeries,1);
                    nCoefficients <- nCoefficients + 1;
                }
                else if(dampedType=="i"){
                    dampedValue <- matrix(B[nCoefficients+(1:nSeries)],nSeries,1);
                    nCoefficients <- nCoefficients + nSeries;
                }
            }

            ### Transition matrix
            if(any(transitionType==c("i","d","c")) & damped){
                for(i in 1:nSeries){
                    matF[c(1:nComponentsNonSeasonal)+nComponentsNonSeasonal*(i-1),
                         nComponentsNonSeasonal+nComponentsNonSeasonal*(i-1)] <- dampedValue[i];
                }
            }
            if(transitionType=="d"){
                # Fill in the other values of F with some values
                nCoefficientsBuffer <- (nSeries-1)*nComponentsNonSeasonal^2;

                for(i in 1:nSeries){
                    matF[c(1:nComponentsNonSeasonal)+nComponentsNonSeasonal*(i-1),
                         setdiff(c(1:(nSeries*nComponentsNonSeasonal)),
                                 c(1:nComponentsNonSeasonal)+nComponentsNonSeasonal*(i-1))
                         ] <- B[nCoefficients+c(1:nCoefficientsBuffer)];
                    nCoefficients <- nCoefficients + nCoefficientsBuffer;
                }
            }

            ### Measurement matrix
            # Needs to be filled in with dampedValue even if dampedValue has been provided by a user
            if(damped){
                for(i in 1:nSeries){
                    matW[i,nComponentsNonSeasonal+nComponentsNonSeasonal*(i-1)] <- dampedValue[i];
                }
            }

            ### Vector of states
            ## Deal with non-seasonal part of the vector of states
            if(initialEstimate){
                initialPlaces <- nComponentsNonSeasonal*(c(1:nSeries)-1)+1;
                if(Ttype!="N"){
                    initialPlaces <- c(initialPlaces,nComponentsNonSeasonal*(c(1:nSeries)-1)+2);
                    initialPlaces <- sort(initialPlaces);
                }
                if(initialType=="i"){
                    initialValue <- matrix(B[nCoefficients+c(1:(nComponentsNonSeasonal*nSeries))],
                                           nComponentsNonSeasonal * nSeries,1);
                    nCoefficients <- nCoefficients + nComponentsNonSeasonal*nSeries;
                }
                else if(initialType=="c"){
                    initialValue <- matrix(B[nCoefficients+c(1:nComponentsNonSeasonal)],nComponentsNonSeasonal * nSeries,1);
                    nCoefficients <- nCoefficients + nComponentsNonSeasonal;
                }
                matVt[initialPlaces,1:lagsModelMax] <- rep(initialValue,lagsModelMax);
            }

            ## Deal with seasonal part of the vector of states
            if(modelIsSeasonal & initialSeasonEstimate){
                matVt[nComponentsNonSeasonal*nSeries+1,1:lagsModelMax] <- matrix(B[nCoefficients+c(1:lagsModelMax)],1,lagsModelMax,byrow=TRUE);
                nCoefficients <- nCoefficients + lagsModelMax;
            }
        }

        return(list(matVt=matVt,matF=matF,matG=matG,matW=matW,dampedValue=dampedValue));
    }

    ##### Basic estimation function for ves() #####
    estimatorVES <- function(...){
        environment(creatorVES) <- environment();
        environment(initialiserVES) <- environment();
        environment(logLikVES) <- environment();
        environment(CF) <- environment();
        environment(fillerVES) <- environment();
        list2env(architectorVES(Etype, Ttype, Stype, damped, nSeries, lags),environment());
        list2env(creatorVES(),environment());

        if(is.null(B) && is.null(ub) && is.null(lb)){
            BList <- initialiserVES(Ttype,Stype,lagsModelMax,nComponentsAll,nComponentsNonSeasonal,nSeries);
            B <- BList$B;

            if(any((B>=BList$BUpper),(B<=BList$BLower))){
                B[B>=BList$BUpper] <- BList$BUpper[B>=BList$BUpper] - 0.01;
                B[B<=BList$BLower] <- BList$BLower[B<=BList$BLower] + 0.01;
            }
        }
        else{
            BList <- initialiserVES(Ttype,Stype,lagsModelMax,nComponentsAll,nComponentsNonSeasonal,nSeries);
            if(is.null(B)){
                B <- BList$B;
            }
            if(!is.null(lb)){
                BList$Blower <- lb;
            }
            if(!is.null(ub)){
                BList$BUpper <- ub;
            }
        }

        maxevalUsed <- maxeval;
        if(is.null(maxeval)){
            maxevalUsed <- length(B) * 40;
        }

        print_level_hidden <- print_level;
        if(print_level==41){
            print_level[] <- 0;
        }

        lossNew <- loss;
        # Change loss to "GV" if it is an additive case (to speed things up)
        # In this case minimum of GV coincides with the MLE of MVNorm
        if(loss=="likelihood" && Etype=="A"){
            lossNew[] <- "GV";
        }

        # Parameters are chosen to speed up the optimisation process and have decent accuracy
        res <- nloptr(B, CF, lb=BList$BLower, ub=BList$BUpper,
                      opts=list(algorithm=algorithm1, xtol_rel=xtol_rel1, maxeval=maxevalUsed, print_level=print_level),
                      loss=lossNew, Etype=Etype, Ttype=Ttype, damped=damped,
                      nComponentsNonSeasonal=nComponentsNonSeasonal, nComponentsAll=nComponentsAll, lagsModelMax=lagsModelMax);
        B <- res$solution;

        # This is just in case something went out of the bounds
        if(any((B>=BList$BUpper),(B<=BList$BLower))){
            BList$BUpper[B>=BList$BUpper] <- B[B>=BList$BUpper] + 1;
            BList$BLower[B<=BList$BLower] <- B[B<=BList$BLower] - 1;
        }

        if(print_level_hidden>0){
            print(res);
        }

        res2 <- nloptr(B, CF, lb=BList$BLower, ub=BList$BUpper,
                       opts=list(algorithm=algorithm2, xtol_rel=xtol_rel2, maxeval=maxevalUsed, print_level=print_level),
                       loss=lossNew, Etype=Etype, Ttype=Ttype, damped=damped,
                       nComponentsNonSeasonal=nComponentsNonSeasonal, nComponentsAll=nComponentsAll, lagsModelMax=lagsModelMax);
        # This condition is needed in order to make sure that we did not make the solution worse
        if(res2$objective <= res$objective){
            res <- res2;
        }
        B <- res$solution;

        if(print_level_hidden>0){
            print(res);
        }

        if(all(B==BList$B) & modelDo=="estimate"){
            if(persistenceEstimate){
                warning(paste0("Failed to optimise the model ETS(", modelCurrent,
                               "). Try different initialisation maybe?\nAnd check all the messages and warnings...",
                               "If you did your best, but the optimiser still fails, report this to the maintainer, please."),
                        call.=FALSE);
            }
        }
        names(B) <- BList$BNames;

        # This is needed for AICc / BICc
        nParamPerSeries <- length(B) / nSeries;

        # First part is for the covariance matrix
        if(loss=="likelihood"){
            nParam <- nSeries * (nSeries + 1) / 2 + length(B);
        }
        else{
            nParam <- nSeries + length(B);
        }

        # likelihood and ICs
        logLikVESValue <- structure(logLikVES(B=B,loss=loss,Etype=Etype),
                                    nobs=obsInSample,df=nParam,class="logLik");
        ICs <- ICsVES(B=B, logLikVESValue=logLikVESValue, nSeries=nSeries,
                      nParamPerSeries=nParamPerSeries, obsInSample=obsInSample);

        # If this is a special case, recalculate CF to get the proper loss value
        if(loss=="likelihood" && Etype=="A"){
            res$objective <- CF(B, loss=loss, Etype=Etype, Ttype=Ttype, damped=damped,
                                nComponentsNonSeasonal=nComponentsNonSeasonal, nComponentsAll=nComponentsAll,
                                lagsModelMax=lagsModelMax);
        }

        # Write down Fisher Information if needed
        if(FI){
            environment(vLikelihoodFunction) <- environment();
            FI <- -numDeriv::hessian(vLikelihoodFunction,B);
            rownames(FI) <- BList$BNames;
            colnames(FI) <- BList$BNames;
        }

        return(list(ICs=ICs,objective=res$objective,B=B,nParam=nParam,logLik=logLikVESValue,FI=FI));
    }

    ##### Function selects ETS components #####
    selectorVES <- function(silent=TRUE,...){
        environment(estimatorVES) <- environment();

        # If the pool is not provided, form it
        if(is.null(modelsPool)){
            if(!silent){
                cat("Estimating models... ");
            }

            # Define the whole pool
            if(!allowMultiplicative){
                poolErrors <- c("A");
                poolTrends <- c("N","A","Ad");
                poolSeasonals <- c("N","A");
            }
            else{
                poolErrors <- c("A","M");
                poolTrends <- c("N","A","Ad","M","Md");
                poolSeasonals <- c("N","A","M");
            }

            # If Etype is not P, then check on additive errors
            if(Etype!="P"){
                poolErrors <- Etype;
            }

            # If Ttype is not P, then create a pool with specified type
            if(Ttype!="P"){
                if(Ttype=="X"){
                    poolTrends <- c("N","A","Ad");
                }
                else if(Ttype=="Y"){
                    poolTrends <- c("N","M","Md");
                }
                else{
                    if(damped){
                        poolTrends <- paste0(Ttype,"d");
                    }
                    else{
                        poolTrends <- Ttype;
                    }
                }
            }

            # If Stype is not P, then create specific pools
            if(Stype!="P"){
                if(Stype=="X"){
                    poolSeasonals <- c("N","A");
                }
                else if(Stype=="Y"){
                    poolSeasonals <- c("N","M");
                }
                else{
                    poolSeasonals <- Stype;
                }
            }

            # Create pool of models
            modelsPool <- paste0(rep(poolErrors,each=length(poolTrends)*length(poolSeasonals)),
                                 poolTrends,
                                 rep(poolSeasonals,each=length(poolTrends)));
        }

        # Leave only pure models
        modelsPool <- modelsPool[modelsPool %in% c("ANN","AAN","AAdN","ANA","AAA","AAdA",
                                                   "MNN","MMN","MMdN","MNM","MMM","MMdM")];

        modelsNumber <- length(modelsPool);
        vesModels <- vector("list", length(modelsPool));
        for(i in 1:modelsNumber){
            if(!silent){
                if(i>1){
                    cat(paste0(rep("\b",nchar(round((i-1)/modelsNumber,2)*100)+1),collapse=""));
                }
                cat(round(i/modelsNumber,2)*100,"\b%");
            }

            modelCurrent <- modelsPool[i];
            Etype <- substring(modelCurrent,1,1);
            Ttype <- substring(modelCurrent,2,2);
            if(nchar(modelCurrent)==4){
                damped <- TRUE;
                Stype <- substring(modelCurrent,4,4);
            }
            else{
                damped <- FALSE;
                Stype <- substring(modelCurrent,3,3);
            }
            vesModels[[i]] <- estimatorVES(ParentEnvironment=environment());
        }
        if(!silent){
            cat("\n");
        }

        # Prepare the return of the best model
        vesModelsICs <- sapply(vesModels,"[[","ICs");
        colnames(vesModelsICs) <- modelsPool;
        iBest <- which.min(vesModelsICs[ic,]);
        vesModels[[iBest]]$model <- modelsPool[iBest];
        vesModels[[iBest]]$Etype <- substring(modelsPool[iBest],1,1);
        vesModels[[iBest]]$Ttype <- substring(modelsPool[iBest],2,2);
        if(nchar(modelsPool[iBest])==4){
            vesModels[[iBest]]$damped <- TRUE;
            vesModels[[iBest]]$Stype <- substring(modelsPool[iBest],4,4);
        }
        else{
            vesModels[[iBest]]$damped <- FALSE;
            vesModels[[iBest]]$Stype <- substring(modelsPool[iBest],3,3);
        }
        vesModels[[iBest]]$ICsAll <- vesModelsICs;
        vesModels[[iBest]]$ICs <- vesModelsICs[,iBest];
        # Rename "objective" into "cfObjective"
        names(vesModels[[iBest]])[names(vesModels[[iBest]])=="objective"] <- "cfObjective";
        return(vesModels[[iBest]]);
    }


    ##### Function constructs the VES function #####
    callerVES <- function(silent=FALSE,...){
        if(modelDo=="estimate"){
            environment(estimatorVES) <- environment();
            res <- estimatorVES(ParentEnvironment=environment());
            listToReturn <- list(Etype=Etype,Ttype=Ttype,Stype=Stype,damped=damped,
                                 cfObjective=res$objective,B=res$B,ICs=res$ICs,
                                 ICsAll=res$ICs,nParam=res$nParam,logLik=res$logLik,FI=res$FI);

            return(listToReturn);
        }
        else if(modelDo=="select"){
            return(selectorVES(silent=silent,ParentEnvironment=environment()));
        }
        else{
            environment(CF) <- environment();
            environment(creatorVES) <- environment();
            environment(logLikVES) <- environment();
            list2env(architectorVES(Etype, Ttype, Stype, damped, nSeries, lags),environment());
            list2env(creatorVES(),environment());

            B <- c(persistenceValue);
            BNames <- paste0("Persistence",c(1:length(persistenceValue)));
            if(damped){
                B <- c(B,dampedValue);
                BNames <- c(BNames,paste0("phi",c(1:length(dampedValue))));
            }
            if(transitionType=="d"){
                transitionLength <- length(B);
                # Write values from the rest of transition matrix
                for(i in 1:nSeries){
                    B <- c(B, c(transitionValue[c(1:nComponentsAll)+nComponentsAll*(i-1),
                                                setdiff(c(1:nSeries*nComponentsAll),
                                                        c(1:nComponentsAll)+nComponentsAll*(i-1))]));
                }
                transitionLength <- length(B) - transitionLength;
                BNames <- c(BNames,paste0("transition",c(1:transitionLength)));
            }
            B <- c(B,initialValue);
            BNames <- c(BNames,paste0("initial",c(1:length(initialValue))));
            if(Stype!="N"){
                B <- c(B,initialSeasonValue);
                BNames <- c(BNames,paste0("initialSeason",c(1:length(initialSeasonValue))));
            }
            names(B) <- BNames;

            cfObjective <- CF(B, loss=loss, Etype=Etype, Ttype=Ttype, damped=damped,
                              nComponentsNonSeasonal=nComponentsNonSeasonal, nComponentsAll=nComponentsAll,
                              lagsModelMax=lagsModelMax);

            # Number of parameters
            # First part is for the covariance matrix
            if(loss=="likelihood"){
                nParam <- nSeries * (nSeries + 1) / 2;
            }
            else if(loss=="diagonal"){
                nParam <- nSeries;
            }
            else{
                nParam <- nSeries;
            }

            # likelihood and ICs
            logLikVESValue <- structure(logLikVES(B=B,loss=loss,Etype=Etype),
                                        nobs=obsInSample,df=nParam,class="logLik");
            ICs <- setNames(c(AIC(logLikVESValue), AICc(logLikVESValue), BIC(logLikVESValue), BICc(logLikVESValue)),
                            c("AIC","AICc","BIC","BICc"));

            # Write down Fisher Information if needed
            if(FI){
                environment(vLikelihoodFunction) <- environment();
                FI <- -numDeriv::hessian(vLikelihoodFunction,B);
                rownames(FI) <- BNames;
                colnames(FI) <- BNames;
            }

            listToReturn <- list(Etype=Etype,Ttype=Ttype,Stype=Stype,damped=damped,
                                 cfObjective=cfObjective,B=B,ICs=ICs,
                                 ICsAll=ICs,nParam=nParam,logLik=logLikVESValue,FI=FI);
            return(listToReturn);
        }
    }

    ##### Preset yFitted, yForecast, errors and basic parameters #####
    yFitted <- matrix(NA,nSeries,obsInSample);
    yForecast <- matrix(NA,nSeries,h);
    errors <- matrix(NA,nSeries,obsInSample);
    rownames(yFitted) <- rownames(yForecast) <- rownames(errors) <- dataNames;

    ##### Define modelDo #####
    if(!any(persistenceEstimate, transitionEstimate, dampedEstimate, initialEstimate, initialSeasonEstimate)){
        modelDo <- "nothing";
        modelCurrent <- model;
        bounds <- "n";
    }

    ##### Now do estimation and model selection #####
    environment(creatorVES) <- environment();
    environment(fillerVES) <- environment();
    environment(vssFitter) <- environment();
    environment(vssForecaster) <- environment();

    ##### Fit the model and produce forecast #####
    list2env(callerVES(silent=silent),environment());
    list2env(architectorVES(Etype, Ttype, Stype, damped, nSeries, lags),environment());
    list2env(creatorVES(),environment());
    list2env(fillerVES(matVt,matF,matG,matW,B,Ttype,damped,
                       nComponentsNonSeasonal,nComponentsAll,lagsModelMax),environment());

    if(damped){
        model <- paste0(Etype,Ttype,"d",Stype);
    }
    else{
        model <- paste0(Etype,Ttype,Stype);
    }

    statesNames <- rownames(matVt);
    vssFitter(ParentEnvironment=environment());
    rownames(matVt) <- statesNames;
    vssForecaster(ParentEnvironment=environment());

    ##### Write down persistence, transition, initials etc #####
    # Write down the persistenceValue, transitionValue, initialValue, initialSeasonValue

    persistenceNames <- "level";
    if(Ttype!="N"){
        persistenceNames <- c(persistenceNames,"trend");
    }
    if(Stype!="N"){
        persistenceNames <- c(persistenceNames,"seasonal");
    }
    if(persistenceEstimate){
        persistenceValue <- matG;
        if(persistenceType=="c"){
            parametersNumber[1,1] <- parametersNumber[1,1] + nComponentsAll;
        }
        else if(persistenceType=="i"){
            if(seasonalType=="i"){
                parametersNumber[1,1] <- parametersNumber[1,1] + nSeries*nComponentsAll;
            }
            else{
                parametersNumber[1,1] <- parametersNumber[1,1] + nSeries*nComponentsNonSeasonal+nSeries;
            }
        }
        else if(persistenceType=="s"){
            if(seasonalType=="i"){
                parametersNumber[1,1] <- parametersNumber[1,1] + nSeries*(nComponentsAll-1)+1;
            }
            else{
                parametersNumber[1,1] <- parametersNumber[1,1] + nSeries*nComponentsNonSeasonal+1;
            }
        }
        else{
            parametersNumber[1,1] <- parametersNumber[1,1] + length(matG);
        }
    }
    if(seasonalType=="i"){
        rownames(persistenceValue) <- paste0(rep(dataNames,each=nComponentsAll), "_", persistenceNames);
    }
    else{
        rownames(persistenceValue) <- c(paste0(rep(dataNames,each=nComponentsNonSeasonal), "_",
                                               persistenceNames[-nComponentsAll]),
                                        persistenceNames[nComponentsAll]);
    }
    colnames(persistenceValue) <- dataNames;

    # This is needed anyway for the reusability of the model
    transitionValue <- matF;
    if(transitionEstimate){
        if(seasonalType=="i"){
            parametersNumber[1,1] <- parametersNumber[1,1] + (nSeries-1)*nSeries*nComponentsAll^2;
        }
        else{
            parametersNumber[1,1] <- parametersNumber[1,1] + (nSeries-1)*nSeries*nComponentsNonSeasonal^2;
        }
    }
    colnames(transitionValue) <- rownames(persistenceValue);
    rownames(transitionValue) <- rownames(persistenceValue);

    if(damped){
        rownames(dampedValue) <- dataNames;
        if(dampedEstimate){
            parametersNumber[1,1] <- parametersNumber[1,1] + length(unique(as.vector(dampedValue)));
        }
    }

    rownames(matW) <- dataNames;
    colnames(matW) <- rownames(persistenceValue);

    if(seasonalType=="i"){
        initialPlaces <- nComponentsAll*(c(1:nSeries)-1)+1;
        initialNames <- "level";
        if(Ttype!="N"){
            initialPlaces <- c(initialPlaces,nComponentsAll*(c(1:nSeries)-1)+2);
            initialPlaces <- sort(initialPlaces);
            initialNames <- c(initialNames,"trend");
        }
        if(initialEstimate){
            initialValue <- matrix(matVt[initialPlaces,lagsModelMax],nComponentsNonSeasonal*nSeries,1);
            parametersNumber[1,1] <- parametersNumber[1,1] + length(unique(as.vector(initialValue)));
        }
    }
    else{
        initialNames <- "level";
        if(Ttype!="N"){
            initialNames <- c(initialNames,"trend");
        }
        if(initialEstimate){
            initialValue <- matrix(matVt[1:(nComponentsNonSeasonal*nSeries),lagsModelMax],
                                   nComponentsNonSeasonal*nSeries,1);
            parametersNumber[1,1] <- parametersNumber[1,1] + length(unique(as.vector(initialValue)));
        }
    }
    rownames(initialValue) <- paste0(rep(dataNames,each=nComponentsNonSeasonal), "_", initialNames);

    if(modelIsSeasonal){
        if(seasonalType=="i"){
            if(initialSeasonEstimate){
                initialPlaces <- nComponentsAll*(c(1:nSeries)-1)+nComponentsAll;
                initialSeasonValue <- matrix(matVt[initialPlaces,1:lagsModelMax],nSeries,lagsModelMax);
                parametersNumber[1,1] <- parametersNumber[1,1] + length(unique(as.vector(initialSeasonValue)));
            }
            rownames(initialSeasonValue) <- dataNames;
        }
        else{
            initialSeasonValue <- matrix(matVt[nComponentsNonSeasonal*nSeries+1,1:lagsModelMax],1,lagsModelMax);
            parametersNumber[1,1] <- parametersNumber[1,1] + lagsModelMax;
            rownames(initialSeasonValue) <- "Common";
        }
        colnames(initialSeasonValue) <- paste0("Seasonal",c(1:lagsModelMax));
    }

    if(!is.matrix(yForecast)){
        yForecast <- as.matrix(yForecast,nSeries,h);
    }
    if(any(yClasses=="ts")){
        yFitted <- ts(t(yFitted), start=yStart, frequency=yFrequency);
        yForecast <- ts(t(yForecast), start=yForecastStart, frequency=yFrequency);
    }
    else{
        yFitted <- zoo(t(yFitted), order.by=yInSampleIndex);
        yForecast <- zoo(t(yForecast), order.by=yForecastIndex);
    }
    colnames(yForecast) <- dataNames;

    if(loss=="likelihood"){
        parametersNumber[1,1] <- parametersNumber[1,1] + nSeries * (nSeries + 1) / 2;
    }
    else if(loss=="diagonal"){
        parametersNumber[1,1] <- parametersNumber[1,1] + nSeries;
    }
    else{
        parametersNumber[1,1] <- parametersNumber[1,1] + nSeries;
    }

    parametersNumber[1,4] <- sum(parametersNumber[1,1:3]);
    parametersNumber[2,4] <- sum(parametersNumber[2,1:3]);

    ##### Now let's deal with the holdout #####
    if(holdout){
        if(any(yClasses=="ts")){
            yHoldout <- ts(data[(obsInSample+1):obsAll,,drop=FALSE], start=yForecastStart, frequency=yFrequency);
        }
        else{
            yHoldout <- zoo(data[(obsInSample+1):obsAll,,drop=FALSE], order.by=yForecastIndex);
        }
        colnames(yHoldout) <- dataNames;

        measureFirst <- measures(yHoldout[,1],yForecast[,1],yInSample[1,]);
        errorMeasures <- matrix(NA,nSeries,length(measureFirst));
        rownames(errorMeasures) <- dataNames;
        colnames(errorMeasures) <- names(measureFirst);
        errorMeasures[1,] <- measureFirst;
        for(i in 2:nSeries){
            errorMeasures[i,] <- measures(yHoldout[,i],yForecast[,i],yInSample[i,]);
        }
    }
    else{
        yHoldout <- NULL;
        errorMeasures <- NA;
    }

    modelname <- "VES";
    modelname <- paste0(modelname,"(",model,")");

    if(occurrence!="n"){
        modelname <- paste0("i",modelname);
    }

    if(modelIsSeasonal){
        submodelName <- "[";
        if(seasonalType=="c"){
            submodelName[] <- paste0(submodelName,"C");
        }
        else{
            submodelName[] <- paste0(submodelName,"I");
        }

        if(persistenceType=="i"){
            submodelName[] <- paste0(submodelName,"I");
        }
        else if(persistenceType=="c"){
            submodelName[] <- paste0(submodelName,"CA");
        }
        else if(persistenceType=="s"){
            submodelName[] <- paste0(submodelName,"CS");
        }
        else if(persistenceType=="d"){
            submodelName[] <- paste0(submodelName,"D");
        }

        if(initialSeasonType=="i"){
            submodelName[] <- paste0(submodelName,"I");
        }
        else{
            submodelName[] <- paste0(submodelName,"C");
        }
        submodelName[] <- paste0(submodelName,"]");
        modelname[] <- paste0(modelname,submodelName);
    }

    ##### Print output #####
    if(!silent){
        if(any(abs(eigen(matF - matG %*% matW)$values)>(1 + 1E-10))){
            warning(paste0("Model VES(",model,") is unstable! ",
                           "Use a different value of 'bounds' parameter to address this issue!"),
                    call.=FALSE);
        }
    }

    ##### Make a plot #####
    # This is a temporary solution
    if(!silent){
        pages <- ceiling(nSeries / 5);
        perPage <- ceiling(nSeries / pages);
        packs <- c(seq(1, nSeries+1, perPage));
        if(packs[length(packs)]<nSeries+1){
            packs <- c(packs,nSeries+1);
        }
        parDefault <- par(no.readonly=TRUE);
        on.exit(par(parDefault));
        for(j in 1:pages){
            par(mar=c(4,4,2,1),mfcol=c(perPage,1));
            for(i in packs[j]:(packs[j+1]-1)){
                # if(any(intervalType==c("u","i","l"))){
                #     plotRange <- range(min(data[,i],yForecast[,i],yFitted[,i],PI[,i*2-1]),
                #                        max(data[,i],yForecast[,i],yFitted[,i],PI[,i*2]));
                # }
                # else{
                    plotRange <- range(min(data[,i],yForecast[,i],yFitted[,i]),
                                       max(data[,i],yForecast[,i],yFitted[,i]));
                # }
                plot(data[,i],main=paste0(modelname," on ",dataNames[i]),ylab="Y",
                     ylim=plotRange, xlim=range(time(data[,i])[1],time(yForecast)[max(h,1)]),
                     type="l");
                lines(yFitted[,i],col="purple",lwd=2,lty=2);
                if(h>1){
                    # if(any(intervalType==c("u","i","l"))){
                    #     lines(PI[,i*2-1],col="darkgrey",lwd=3,lty=2);
                    #     lines(PI[,i*2],col="darkgrey",lwd=3,lty=2);
                    #
                    #     polygon(c(seq(yDeltat*(yForecastStart[2]-1)+yForecastStart[1],
                    #                   yDeltat*(end(yForecast)[2]-1)+end(yForecast)[1],yDeltat),
                    #               rev(seq(yDeltat*(yForecastStart[2]-1)+yForecastStart[1],
                    #                       yDeltat*(end(yForecast)[2]-1)+end(yForecast)[1],yDeltat))),
                    #             c(as.vector(PI[,i*2]), rev(as.vector(PI[,i*2-1]))), col="lightgray",
                    #             border=NA, density=10);
                    # }
                    lines(yForecast[,i],col="blue",lwd=2);
                }
                else{
                    # if(any(intervalType==c("u","i","l"))){
                    #     points(PI[,i*2-1],col="darkgrey",lwd=3,pch=4);
                    #     points(PI[,i*2],col="darkgrey",lwd=3,pch=4);
                    # }
                    points(yForecast[,i],col="blue",lwd=2,pch=4);
                }
                abline(v=yForecastStart-yDeltat,col="red",lwd=2);
            }
        }
    }

    ##### Return values #####
    model <- list(model=modelname,timeElapsed=Sys.time()-startTime,
                  states=NA,persistence=persistenceValue,transition=transitionValue,
                  measurement=matW, phi=dampedValue, B=B,
                  lagsAll=lagsModel,
                  initialType=initialType,initial=initialValue,initialSeason=initialSeasonValue,
                  nParam=parametersNumber, imodel=ovesModel,
                  data=NA,fitted=yFitted,holdout=yHoldout,residuals=NA,Sigma=Sigma,
                  forecast=yForecast,
                  # PI=PI,interval=intervalType,level=level,
                  ICs=ICs,ICsAll=ICsAll,logLik=logLik,
                  lossValue=cfObjective,loss=loss,lossFunction=lossFunction,
                  accuracy=errorMeasures,FI=FI);
    # Produce proper objects and return them
    if(any(yClasses=="ts")){
        model$states <- ts(t(matVt), start=(time(data)[1] - yDeltat*lagsModelMax), frequency=yFrequency)
        model$data <- ts(t(yInSample), start=yStart, frequency=yFrequency);
        model$residuals <- ts(t(errors), start=yStart, frequency=yFrequency);
    }
    else{
        yStatesIndex <- yInSampleIndex[1] - lagsModelMax*diff(tail(yInSampleIndex,2)) +
            c(1:lagsModelMax-1)*diff(tail(yInSampleIndex,2));
        yStatesIndex <- c(yStatesIndex, yInSampleIndex);
        model$states <- zoo(t(matVt), order.by=yStatesIndex);
        model$data <- zoo(t(yInSample), order.by=yInSampleIndex);
        model$residuals <- zoo(t(errors), order.by=yInSampleIndex);
    }

    return(structure(model,class=c("legion","smooth")));
}
