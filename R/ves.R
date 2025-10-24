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
                         "transitionValue","dampedValue",
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
#' Pure multiplicative models are done as additive ones applied to log(y).
#'
#' Also \code{model} can accept a previously estimated VES model and use all its
#' parameters.
#' @param lags The lags of the model. Needed for seasonal models.
#' @param phi If vector or a value is provided here, then it is used by the
#' model.
#' @param initial Can be either character or a vector / matrix of initial states.
#' If it is character, then it can be \code{"optimal"}, or \code{"backcasting"}. The
#' former means that the initial states will be estimated in the optimisation, while
#' the latter means that they will be produced via backasting.
#' @param persistence Persistence matrix \eqn{G}, containing smoothing
#' parameters. Can be:
#' \itemize{
#' \item \code{"independent"} - each series has its own smoothing parameters
#' and no interactions are modelled (all the other values in the matrix are set
#' to zero);
#' \item \code{"dependent"} - each series has its own smoothing parameters, and
#' interactions between the series are modelled (the whole matrix is estimated);
#' \item provided by user as a vecor or a matrix.
#' }
#' @param transition Transition matrix \eqn{F}. Can be:
#' \itemize{
#' \item \code{"independent"} - each series has its own preset transition matrix
#' and no interactions are modelled (all the other values in the matrix are set
#' to zero);
#' \item \code{"dependent"} - each series has its own transition matrix, and
#' interactions between the series are modelled (the whole matrix is estimated). The
#' estimated model behaves similar to VAR in this case;
#' \item provided by user as a vector or as a matrix. The value is used by the model.
#' }
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
                persistence=c("independent","dependent"),
                transition=c("independent","dependent"), phi=NULL,
                initial=c("backcasting","optimal"),
                loss=c("likelihood","diagonal","trace","log-trace"),
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
        initialType <- model$initialType;
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
        if(bounds=="admissible"){
            # eigenValues <- eigen(elements$matF - elements$matG %*% elements$matW, only.values=TRUE, symmetric=FALSE)$values;
            # if(nComponentsAll>10){
            #     eigenValues <- discounter(elements$matF, elements$matW, elements$matG, 5);
            # }
            # else{
            #     eigenValues <- eigen(elements$matF - elements$matG %*% elements$matW, only.values=TRUE)$values;
            # }
            eigenValues <- eigenChecker(elements$matF, elements$matW, elements$matG,
                                        nComponentsNonSeasonal, nComponentsAll,
                                        nComponentsSeasonal, modelIsSeasonal, nSeries);
            if(max(abs(eigenValues)>(1 + 1E-50))){
                return(max(abs(eigenValues))*1E+100);
            }
        }

        # Fit the model
        fitting <- vFitterWrap(switch(Etype, "M"=log(yInSample), yInSample),
                               elements$matVt, elements$matF,
                               elements$matW, elements$matG,
                               lagsModel, Etype, Ttype, Stype,
                               ot, initialType=="backcasting",
                               nComponentsTrend);

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
        else if(loss=="log-trace"){
            cfRes <- sum(log(rowSums(fitting$errors^2))) / obsInSample;
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

        ### Persistence matrix
        # Edit KFP: change the parameter space, make all persistence to have mean-reverting behaviour
        # by forcing the off-diagonals to be negative so that the MA coefficient would be positive
        # potentially works for VES(ANN) only for now - dependent persistence matrix
        if(persistenceEstimate){
            if(persistenceType=="independent"){
                persistenceLength <- nComponentsAll*nSeries;
                if(bounds=="usual"){
                    BLower <- c(BLower,rep(0,persistenceLength));
                    BUpper <- c(BUpper,rep(1,persistenceLength));
                }
                else{
                    BLower <- c(BLower,rep(-5,persistenceLength));
                    BUpper <- c(BUpper,rep(5,persistenceLength));
                }
            }
            # Edited by KFP
            else{
                persistenceLength <- nComponentsAll*nSeries^2;
                if(bounds=="usual"){
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
        # Only makes sense when transitionEstimate==FALSE, i.e. independent transition matrix
        if(dampedEstimate){
            B <- c(B,rep(0.95,nSeries));
            BLower <- c(BLower,rep(0,nSeries));
            BUpper <- c(BUpper,rep(1,nSeries));
            BNames <- c(BNames,paste0("phi",c(1:nSeries)));
        }

        ### Transition matrix
        if(transitionEstimate){
            transitionLength <- (nSeries*nComponentsAll^2)*nSeries;
            B <- c(B,rep(0.1,transitionLength));
            BLower <- c(BLower,rep(-1,transitionLength));
            BUpper <- c(BUpper,rep(1,transitionLength));
            BNames <- c(BNames,paste0("transition",c(1:transitionLength)));
        }

        ### Vector of initials
        # Not needed in case of backcasting or provided
        if(initialEstimate && (initialType=="optimal")){
            initialLength <- nComponentsNonSeasonal*nSeries;
            B <- c(B,initialValue);
            BNames <- c(BNames,paste0("initial",c(1:initialLength)));
            BLower <- c(BLower,rep(-Inf,initialLength));
            BUpper <- c(BUpper,rep(Inf,initialLength));

            ### Vector of initial seasonals
            if(modelIsSeasonal){
                initialSeasonLength <- (lagsModelMax-1)*nSeries;
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

        # Here, we use components per series!
        # !!! VETS uses overall number !!!
        nComponentsNonSeasonal <- 1 + (Ttype!="N")*1;
        nComponentsSeasonal <- modelIsSeasonal*sum(lags!=1);
        nComponentsAll <- nComponentsNonSeasonal + nComponentsSeasonal;

        return(list(modelIsTrendy=modelIsTrendy, modelIsSeasonal=modelIsSeasonal, lagsModelMax=lagsModelMax,
                    nComponentsNonSeasonal=nComponentsNonSeasonal, nComponentsSeasonal=nComponentsSeasonal,
                    nComponentsAll=nComponentsAll));
    }

    ##### Basic matrices creator #####
    # This thing returns matVt, matF, matG, matW, dampedValue, initialValue
    # and initialSeasonValue if they are not provided + lagsModel
    creatorVES <- function(...){
        # ellipsis <- list(...);
        # ParentEnvironment <- ellipsis[['ParentEnvironment']];

        #### Persistence matrix ####
        if(persistenceType!="provided"){
            matGBlockAlpha <- diag(nSeries);
            if(modelIsTrendy){
                matGBlockBeta <- diag(nSeries);
            }
            else{
                matGBlockBeta <- NULL;
            }
            if(modelIsSeasonal){
                matGBlockGamma <- diag(nSeries);
            }
            else{
                matGBlockGamma <- NULL;
            }
            matG <- Matrix(rbind(matGBlockAlpha, matGBlockBeta, matGBlockGamma), sparse=TRUE);
        }
        else{
            matG <- Matrix(persistenceValue, sparse=TRUE);
        }


        #### Damping parameter ####
        if(dampedEstimate){
            dampedValue <- 0.95;
        }
        else{
            dampedValue <- 1;
        }


        #### Transition matrix ####
        if(transitionType!="provided"){
            # Preset values in case there's no trend and seasonality
            matFBlock12 <-
            matFBlock21 <-
            matFBlock22 <-
            matFBlock13 <-
            matFBlock31 <-
            matFBlock33 <-
            matFBlock23 <-
            matFBlock32 <- NULL

            if(transitionType=="independent"){
                matFBlock11 <- diag(nSeries)
                # L & T elements
                if(modelIsTrendy){
                    matFBlock12 <- dampedValue*diag(nSeries);
                    matFBlock21 <- matrix(0,nSeries,nSeries);
                    matFBlock22 <- dampedValue*diag(nSeries);
                }

                # L & S elements
                if(modelIsSeasonal){
                    matFBlock13 <- matrix(0,nSeries,nSeries);
                    matFBlock31 <- matrix(0,nSeries,nSeries);
                    matFBlock33 <- diag(nSeries);
                }

                # T & S elements
                if(modelIsTrendy && modelIsSeasonal){
                    matFBlock23 <- matrix(0,nSeries,nSeries);
                    matFBlock32 <- matrix(0,nSeries,nSeries);
                }
            }
            else{
                matFBlock11 <- matrix(0.1, nSeries, nSeries);
                # L & T elements
                if(modelIsTrendy){
                    matFBlock12 <- matFBlock21 <- matFBlock22 <- matFBlock11;
                }

                # L & S elements
                if(modelIsSeasonal){
                    matFBlock13 <- matFBlock31 <- matFBlock33 <- matFBlock11;
                }

                # T & S elements
                if(modelIsTrendy && modelIsSeasonal){
                    matFBlock23 <- matFBlock32 <- matFBlock11;
                }
            }
            matF <- Matrix(rbind(cbind(matFBlock11,matFBlock12,matFBlock13,deparse.level=0),
                                 cbind(matFBlock21,matFBlock22,matFBlock23,deparse.level=0),
                                 cbind(matFBlock31,matFBlock32,matFBlock33,deparse.level=0)), sparse=TRUE);
        }
        else{
            matF <- Matrix(transitionValue, sparse=TRUE);
        }


        #### Blocks for measurement matrix ####
        matWBlock1 <- diag(nSeries);
        if(modelIsTrendy){
            matWBlock2 <- dampedValue*diag(nSeries);
        }
        else{
            matWBlock2 <- NULL;
        }
        if(modelIsSeasonal){
            matWBlock3 <- diag(nSeries);
        }
        else{
            matWBlock3 <- NULL;
        }
        matW <- Matrix(cbind(matWBlock1,matWBlock2,matWBlock3,deparse.level=0), sparse=TRUE);


        #### Vector of states ####
        matVt <- matrix(NA, nComponentsAll * nSeries, obsStates);

        XValues <- rbind(rep(1,obsInSample),c(1:obsInSample));
        initialValue <- switch(Etype, "M"=log(yInSample), yInSample) %*% t(XValues) %*% solve(XValues %*% t(XValues));

        #### !!! This is a temporary fix for log(0) on intermittent demand ####
        if(any(is.nan(initialValue))){
            initialValue[,1] <- rowMeans(yInSample);
            initialValue[,2] <- switch(Ttype, "M"=1, 0);
        }

        j <- 0;
        # Fill in level initials
        # We have nSeries level components
        matVt[j+1:nSeries,1:lagsModelMax] <- initialValueNew01 <- initialValue[,1];
        j[] <- j+nSeries;

        # Fill in trend initials
        if(modelIsTrendy){
        # And nSeries trend components
            matVt[j+1:nSeries,1:lagsModelMax] <- initialValueNew02 <- initialValue[,2];
            j[] <- j+nSeries;
        }
        else{
            initialValueNew02 <- NULL;
        }
        initialValue <- c(initialValueNew01, initialValueNew02);

        if(modelIsSeasonal){
            if(modelIsTrendy && !shortSample){
                # Matrix of linear trend and dummies for seasons
                XValues <- rbind(rep(1,obsInSample),c(1:obsInSample),
                                 matrix(rep(diag(lagsModelMax)[-1,],ceiling(obsInSample/lagsModelMax)),
                                        lagsModelMax-1)[,1:obsInSample]);
                initialSeasonValue <- (switch(Etype, "M"=log(yInSample), yInSample) %*%
                                           t(XValues) %*% solve(XValues %*% t(XValues)))[,-2];
            }
            else{
                # Matrix of dummies for seasons
                XValues <- rbind(rep(1,obsInSample),
                                 matrix(rep(diag(lagsModelMax)[-1,],ceiling(obsInSample/lagsModelMax)),
                                        lagsModelMax-1)[,1:obsInSample]);
                initialSeasonValue <- (switch(Etype, "M"=log(yInSample), yInSample) %*%
                                           t(XValues) %*% solve(XValues %*% t(XValues)));
            }

            # Fix the seasonals based on dummies
            initialSeasonValue[,-1] <- initialSeasonValue[,-1] + initialSeasonValue[,1];
            # Renormalise initials
            initialSeasonValue[] <- initialSeasonValue - rowMeans(initialSeasonValue);

            matVt[j+1:(nComponentsSeasonal*nSeries),1:lagsModelMax] <- initialSeasonValue;

            # Remove one of columns, to preserve degrees of freedom (normalisation)
            initialSeasonValue <- initialSeasonValue[,-1,drop=FALSE];
        }
        else{
            initialSeasonValue <- NULL;
        }

        #### Names of states ####
        statesNames <- paste0(rep("level",nSeries),"_",dataNames);
        if(modelIsTrendy){
            statesNames <- c(statesNames,paste0(rep("trend",nSeries),"_",dataNames));
        }
        if(modelIsSeasonal){
            statesNames <- c(statesNames,paste0(rep("seasonal",nSeries),"_",dataNames));
        }
        # Give proper names to all matrices
        rownames(matVt) <- statesNames;
        rownames(matF) <- colnames(matF) <- statesNames;
        colnames(matW) <- rownames(matG) <- statesNames;
        rownames(matW) <- colnames(matG) <- dataNames;


        ### lagsModel vector
        lagsModel <- matrix(1,nComponentsAll*nSeries,1);
        if(modelIsSeasonal){
            lagsModel[nComponentsNonSeasonal*nSeries + 1:(nComponentsSeasonal*nSeries)] <- lagsModelMax;
        }

        return(list(matVt=matVt, matF=matF, matG=matG, matW=matW, dampedValue=dampedValue,
                    initialValue=initialValue, initialSeasonValue=initialSeasonValue, lagsModel=lagsModel));
    }

    ##### Basic matrices filler #####
    # This thing fills in matVt, matF, matG and matW with values from B and returns the corrected values
    fillerVES <- function(matVt,matF,matG,matW,B,Ttype,damped,
                          nComponentsNonSeasonal,nComponentsAll,lagsModelMax){

        nCoefficients <- 0;
        ### Persistence matrix
        if(persistenceEstimate){
            # Dependent values
            if(persistenceType=="dependent"){
                matG[] <- B[1:(nComponentsAll*nSeries^2)];
                nCoefficients[] <- nComponentsAll*nSeries^2;
            }
            else if(persistenceType=="independent"){
                matG[1:nSeries,1:nSeries] <- B[nCoefficients+1:nSeries]*diag(nSeries);
                nCoefficients[] <- nCoefficients+nSeries;
                if(modelIsTrendy){
                    matG[nCoefficients+1:nSeries,1:nSeries] <- B[nCoefficients+1:nSeries]*diag(nSeries);
                    nCoefficients[] <- nCoefficients+nSeries;
                }
                if(modelIsSeasonal){
                    matG[nCoefficients+1:nSeries,1:nSeries] <- B[nCoefficients+1:nSeries]*diag(nSeries);
                    nCoefficients[] <- nCoefficients+nSeries;
                }
            }
        }

        ### Damping parameter
        if(modelIsTrendy && dampedEstimate && transitionType=="independent"){
            # Damped parameter in the level equations
            matF[1:nSeries,nSeries+1:nSeries] <-
            # Damped parameter in the trend equations
            matF[nSeries+1:nSeries,nSeries+1:nSeries] <-
            # Measurement matrix
            matW[1:nSeries,nSeries+1:nSeries] <- B[nCoefficients+1:nSeries]*diag(nSeries);

            nCoefficients[] <- nCoefficients + nSeries;
        }

        ### Transition matrix
        if(transitionType=="dependent"){
            matF[1:(nComponentsAll*nSeries)^2] <- B[nCoefficients+1:(nComponentsAll*nSeries)^2]
            nCoefficients[] <- nCoefficients + (nComponentsAll*nSeries)^2;
        }

        ### Vector of states
        ## Deal with non-seasonal part of the vector of states
        if(initialEstimate && initialType=="optimal"){
            # Index for the states rows
            k <- 0;

            # Fill in level initials
            matVt[1:nSeries,1:lagsModelMax] <- B[nCoefficients+1:nSeries];
            nCoefficients[] <- nCoefficients + nSeries;
            k[] <- k + nSeries;

            # Fill in trend initials
            if(modelIsTrendy){
                matVt[nSeries+1:nSeries,1:lagsModelMax] <- B[nCoefficients+1:nSeries];
                nCoefficients[] <- nCoefficients + nSeries;
                k[] <- k + nSeries;
            }

            # Fill in seasonals
            if(modelIsSeasonal){
                matVt[k+1:(nComponentsSeasonal*nSeries),2:lagsModelMax] <-
                    matrix(B[nCoefficients+1:(nSeries*(lagsModelMax-1))],
                           nComponentsSeasonal*nSeries,(lagsModelMax-1),byrow=TRUE);

                # Normalise initials
                matVt[k+1:(nComponentsSeasonal*nSeries),1] <-
                    -rowSums(matVt[k+1:(nComponentsSeasonal*nSeries),2:lagsModelMax,drop=FALSE]);
            }
        }

        return(list(matVt=matVt,matF=matF,matG=matG,matW=matW));
    }

    ##### Basic estimation function for ves() #####
    estimatorVES <- function(...){
        environment(creatorVES) <- environment();
        environment(initialiserVES) <- environment();
        environment(logLikVES) <- environment();
        environment(CF) <- environment();
        environment(eigenChecker) <- environment();
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

        if(print_level_hidden>0){
            print(B);
        }

        # If B has zero length, then everything was provided.
        if(length(B>0)){
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
        }
        else{
            res <- list(objective=CF(B,loss=lossNew, Etype=Etype, Ttype=Ttype, damped=damped,
                           nComponentsNonSeasonal=nComponentsNonSeasonal,
                           nComponentsAll=nComponentsAll, lagsModelMax=lagsModelMax),
                        solution=B);
        }

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
    if(!any(persistenceEstimate, transitionEstimate, dampedEstimate, initialEstimate)){
        modelDo <- "nothing";
        modelCurrent <- model;
        bounds <- "none";
    }

    ##### Now do estimation and model selection #####
    environment(creatorVES) <- environment();
    environment(fillerVES) <- environment();
    environment(vssFitter) <- environment();
    environment(vssForecaster) <- environment();

    nComponentsTrend <- nSeries;

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

    parametersNumber[1,1] <- length(B);
#
#     persistenceNames <- "level";
#     if(Ttype!="N"){
#         persistenceNames <- c(persistenceNames,"trend");
#     }
#     if(Stype!="N"){
#         persistenceNames <- c(persistenceNames,"seasonal");
#     }
#
#
#     if(persistenceEstimate){
#         persistenceValue <- matG;
#         if(persistenceType=="c"){
#             parametersNumber[1,1] <- parametersNumber[1,1] + nComponentsAll;
#         }
#         else if(persistenceType=="i"){
#             if(seasonalType=="i"){
#                 parametersNumber[1,1] <- parametersNumber[1,1] + nSeries*nComponentsAll;
#             }
#             else{
#                 parametersNumber[1,1] <- parametersNumber[1,1] + nSeries*nComponentsNonSeasonal+nSeries;
#             }
#         }
#         else if(persistenceType=="s"){
#             if(seasonalType=="i"){
#                 parametersNumber[1,1] <- parametersNumber[1,1] + nSeries*(nComponentsAll-1)+1;
#             }
#             else{
#                 parametersNumber[1,1] <- parametersNumber[1,1] + nSeries*nComponentsNonSeasonal+1;
#             }
#         }
#         else{
#             parametersNumber[1,1] <- parametersNumber[1,1] + length(matG);
#         }
#     }
#     if(seasonalType=="i"){
#         rownames(persistenceValue) <- paste0(rep(dataNames,each=nComponentsAll), "_", persistenceNames);
#     }
#     else{
#         rownames(persistenceValue) <- c(paste0(rep(dataNames,each=nComponentsNonSeasonal), "_",
#                                                persistenceNames[-nComponentsAll]),
#                                         persistenceNames[nComponentsAll]);
#     }
#     colnames(persistenceValue) <- dataNames;
#
#     # This is needed anyway for the reusability of the model
#     transitionValue <- matF;
#     if(transitionEstimate){
#         if(seasonalType=="i"){
#             parametersNumber[1,1] <- parametersNumber[1,1] + (nSeries-1)*nSeries*nComponentsAll^2;
#         }
#         else{
#             parametersNumber[1,1] <- parametersNumber[1,1] + (nSeries-1)*nSeries*nComponentsNonSeasonal^2;
#         }
#     }
#     colnames(transitionValue) <- rownames(persistenceValue);
#     rownames(transitionValue) <- rownames(persistenceValue);
#
#     if(damped){
#         rownames(dampedValue) <- dataNames;
#         if(dampedEstimate){
#             parametersNumber[1,1] <- parametersNumber[1,1] + length(unique(as.vector(dampedValue)));
#         }
#     }
#
#     rownames(matW) <- dataNames;
#     colnames(matW) <- rownames(persistenceValue);
#
#     if(seasonalType=="i"){
#         initialPlaces <- nComponentsAll*(c(1:nSeries)-1)+1;
#         initialNames <- "level";
#         if(Ttype!="N"){
#             initialPlaces <- c(initialPlaces,nComponentsAll*(c(1:nSeries)-1)+2);
#             initialPlaces <- sort(initialPlaces);
#             initialNames <- c(initialNames,"trend");
#         }
#         if(initialEstimate){
#             initialValue <- matrix(matVt[initialPlaces,lagsModelMax],nComponentsNonSeasonal*nSeries,1);
#             parametersNumber[1,1] <- parametersNumber[1,1] + length(unique(as.vector(initialValue)));
#         }
#     }
#     else{
#         initialNames <- "level";
#         if(Ttype!="N"){
#             initialNames <- c(initialNames,"trend");
#         }
#         if(initialEstimate){
#             initialValue <- matrix(matVt[1:(nComponentsNonSeasonal*nSeries),lagsModelMax],
#                                    nComponentsNonSeasonal*nSeries,1);
#             parametersNumber[1,1] <- parametersNumber[1,1] + length(unique(as.vector(initialValue)));
#         }
#     }
#     rownames(initialValue) <- paste0(rep(dataNames,each=nComponentsNonSeasonal), "_", initialNames);
#
#     if(modelIsSeasonal){
#         if(seasonalType=="i"){
#             if(initialSeasonEstimate){
#                 initialPlaces <- nComponentsAll*(c(1:nSeries)-1)+nComponentsAll;
#                 initialSeasonValue <- matrix(matVt[initialPlaces,1:lagsModelMax],nSeries,lagsModelMax);
#                 parametersNumber[1,1] <- parametersNumber[1,1] + length(unique(as.vector(initialSeasonValue)));
#             }
#             rownames(initialSeasonValue) <- dataNames;
#         }
#         else{
#             initialSeasonValue <- matrix(matVt[nComponentsNonSeasonal*nSeries+1,1:lagsModelMax],1,lagsModelMax);
#             parametersNumber[1,1] <- parametersNumber[1,1] + lagsModelMax;
#             rownames(initialSeasonValue) <- "Common";
#         }
#         colnames(initialSeasonValue) <- paste0("Seasonal",c(1:lagsModelMax));
#     }

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

    submodelName <- "-TP["
    if(transitionType=="independent"){
        submodelName[] <- paste0(submodelName,"I,");
    }
    else{
        submodelName[] <- paste0(submodelName,"D,");
    }
    if(persistenceType=="independent"){
        submodelName[] <- paste0(submodelName,"I]");
    }
    else{
        submodelName[] <- paste0(submodelName,"D]");
    }
    modelname[] <- paste0(modelname,submodelName);

    # if(modelIsSeasonal){
    #     submodelName <- "[";
    #     if(seasonalType=="c"){
    #         submodelName[] <- paste0(submodelName,"C");
    #     }
    #     else{
    #         submodelName[] <- paste0(submodelName,"I");
    #     }
    #
    #     if(persistenceType=="i"){
    #         submodelName[] <- paste0(submodelName,"I");
    #     }
    #     else if(persistenceType=="c"){
    #         submodelName[] <- paste0(submodelName,"CA");
    #     }
    #     else if(persistenceType=="s"){
    #         submodelName[] <- paste0(submodelName,"CS");
    #     }
    #     else if(persistenceType=="d"){
    #         submodelName[] <- paste0(submodelName,"D");
    #     }
    #
    #     if(initialSeasonType=="i"){
    #         submodelName[] <- paste0(submodelName,"I");
    #     }
    #     else{
    #         submodelName[] <- paste0(submodelName,"C");
    #     }
    #     submodelName[] <- paste0(submodelName,"]");
    #     modelname[] <- paste0(modelname,submodelName);
    # }

    ##### Print output #####
    eigenValues <- eigenChecker(matF, matW, matG,
                                nComponentsNonSeasonal, nComponentsAll,
                                nComponentsSeasonal, modelIsSeasonal, nSeries)
    if(any(abs(silent)>(1 + 1E-10))){
        warning(paste0("Model VES(",model,") is unstable! ",
                       "Use a different value of 'bounds' parameter to address this issue!"),
                call.=FALSE);
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
                  states=NA,persistence=as.matrix(matG), transition=as.matrix(matF),
                  measurement=as.matrix(matW), phi=dampedValue, B=B,
                  lagsAll=lagsModel,
                  initialType=initialType, initial=initialValue,
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
