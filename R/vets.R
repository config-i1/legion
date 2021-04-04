utils::globalVariables(c("obsInSample","componentsCommonLevel","componentsCommonSeasonal","componentsCommonTrend",
                         "initialValue","initialsCommonLevel","initialsCommonSeasonal","initialsCommonTrend",
                         "modelIsTrendy","nComponentsLevel","nComponentsSeasonal","nComponentsTrend",
                         "nInitialsLevel","nInitialsSeasonal","nInitialsTrend",
                         "nParametersDamped","nParametersLevel","nParametersSeasonal","nParametersTrend",
                         "parametersCommonDamped","parametersCommonLevel","parametersCommonSeasonal","parametersCommonTrend",
                         "allowMultiplicative","modelDo","ICsAll"));

#' Vector ETS-PIC model
#'
#' Function constructs vector ETS model based on VETS-PIC taxonomy and returns
#' forecast, fitted values, errors and matrix of states along with other useful variables.
#'
#' Function estimates vector ETS in the form of the Single Source of Error state space
#' model of the following type:
#'
#' \deqn{
#' \mathbf{y}_{t} = \mathbf{o}_{t} (\mathbf{W} \mathbf{v}_{t-l} + \mathbf{x}_t
#' \mathbf{a}_{t-1} + \mathbf{\epsilon}_{t})
#' }{
#' y_{t} = o_{t} (W v_{t-l} + x_t a_{t-1} + \epsilon_{t})
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
#' Where \eqn{y_{t}} is the vector of time series on observation \eqn{t}, \eqn{o_{t}}
#' is the vector of Bernoulli distributed random variable (in case of normal data it
#' becomes unit vector for all observations), \eqn{\mathbf{v}_{t}} is the matrix of
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
#' The main idea of the function is in imposing restrictions on parameters / initials /
#' components of the model in order to capture the common dynamics between series.
#'
#' In case of multiplicative model, instead of the vector y_t we use its logarithms.
#' As a result the multiplicative model is much easier to work with.
#'
#' For some more information about the model and its implementation, see the
#' vignette: \code{vignette("ves","legion")}
#'
#' @template vssBasicParam
#' @template vssAdvancedParam
#' @template vssIntervals
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
#' Also \code{model} can accept a previously estimated VETS model and use all its
#' parameters.
#' @param lags The lags of the model. Needed for seasonal models.
#' @param parameters The character vector, specifying, which of the parameters
#' should be common between time series. This includes smoothing parameters for
#' \code{"level"}, \code{"trend"}, \code{"seasonal"} components and \code{"damped"}
#' trend parameter. If \code{parameters="none"}, then all parameters are set to be
#' individual. An example is the model with all parameters being common:
#' \code{parameters=c("level","trend","seasonal","damped")}. The order is not important
#' and the first letters can be used instead of the full words as well.
#' @param initials The character vector, specifying, which of the initial values of
#' components should be common. This can be \code{"level"}, \code{"trend"} and / or
#' \code{"seasonal"}, setting initials of respective components to be common. This
#' can also be \code{"none"}, making the initials individual for all series. An
#' example is the model with only seasonal initials bein common:
#' \code{initials="seasonal"}. The order is not important, and the first letters can
#' be used instead of the full words.
#' @param components The character vector, specifying, which of the components
#' components should be shared between time series. This can be \code{"level"},
#' \code{"trend"} and / or \code{"seasonal"}, setting respective components to be
#' shared. This can also be \code{"none"}, making them individual for all series.
#' The order is not important, and the first letters can be used instead of the full
#' words. Please, note that making components common automatically sets the
#' respective \code{initials} common as well.
#' @param ... Other non-documented parameters. For example \code{FI=TRUE} will
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
#' \item \code{B} - The vector of all the estimated coefficients;
#' \item \code{initial} - The initial values of the non-seasonal components;
#' \item \code{initialSeason} - The initial values of the seasonal components;
#' \item \code{nParam} - The number of estimated parameters;
#' \item \code{occurrence} - The occurrence model estimated with VETS;
#' \item \code{y} - The matrix with the original data;
#' \item \code{fitted} - The matrix of the fitted values;
#' \item \code{holdout} - The matrix with the holdout values (if \code{holdout=TRUE} in
#' the estimation);
#' \item \code{residuals} - The matrix of the residuals of the model;
#' \item \code{Sigma} - The covariance matrix of the errors (estimated with the correction
#' for the number of degrees of freedom);
#' \item \code{forecast} - The matrix of point forecasts;
#' \item \code{PI} - The bounds of the prediction interval;
#' \item \code{interval} - The type of the constructed prediction interval;
#' \item \code{level} - The level of the confidence for the prediction interval;
#' \item \code{ICs} - The values of the information criteria;
#' \item \code{logLik} - The log-likelihood function;
#' \item \code{lossValue} - The value of the loss function;
#' \item \code{loss} - The type of the used loss function;
#' \item \code{accuracy} - the values of the error measures. Currently not available.
#' \item \code{FI} - Fisher information if user asked for it using \code{FI=TRUE}.
#' }
#' @seealso \code{\link[legion]{ves}, \link[smooth]{es}, \link[smooth]{adam}}
#'
#' @examples
#'
#' Y <- ts(cbind(rnorm(100,100,10),rnorm(100,75,8)),frequency=12)
#'
#' # The simplest model applied to the data with the default values
#' vets(Y,model="ANN",h=10,holdout=TRUE)
#'
#' # Multiplicative damped trend model with common parameters
#' # and initial seasonal indices
#' vets(Y,model="MMdM",h=10,holdout=TRUE,parameters=c("l","t","s","d"),
#'      initials="seasonal")
#'
#' # Automatic selection of ETS components
#' vets(Y, model="PPP", h=10, holdout=TRUE, initials="seasonal")
#'
#' @export
vets <- function(y, model="ANN", lags=c(frequency(y)),
                 parameters=c("level","trend","seasonal","damped"),
                 initials=c("seasonal"), components=c("none"),
                 loss=c("likelihood","diagonal","trace"),
                 ic=c("AICc","AIC","BIC","BICc"), h=10, holdout=FALSE,
                 interval=c("none","conditional","unconditional","individual","likelihood"), level=0.95,
                 occurrence=c("none","fixed","logistic"),
                 bounds=c("admissible","usual","none"),
                 silent=c("all","graph","output","none"), ...){
    # Copyright (C) 2017 - Inf  Ivan Svetunkov

    # Start measuring the time of calculations
    startTime <- Sys.time();

    cumulative <- FALSE;

    # If a previous model provided as a model, write down the variables
    #### This needs to be updated, when VETS works! ####
    if(any(is.legion(model))){
        persistence <- model$persistence;
        transition <- model$transition;
        measurement <- model$measurement;
        initial <- model$initial;
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
    vssInput("vets",ParentEnvironment=environment(),...);

    ##### Basic VETS architector #####
    ### This function will accept Etype, Ttype, Stype and damped and would return:
    # nComponentsNonSeasonal, nComponentsAll, lagsModelMax, modelIsSeasonal, modelIsTrendy, obsStates
    # and all the variables for PIC part of the model
    # This is needed for model selection
    architectorVETS <- function(Etype, Ttype, Stype, damped, nSeries){
        # Binaries for trend and seasonal
        modelIsTrendy <- Ttype!="N";
        modelIsSeasonal <- Stype!="N";

        lagsModelMax <- dataFreq * modelIsSeasonal + 1 * (!modelIsSeasonal);

        # Define the number of rows that should be in the matVt
        obsStates <- max(obsAll + lagsModelMax, obsInSample + 2*lagsModelMax);

        # Common parameters
        parametersCommonLevel <- any(parameters=="l");
        parametersCommonTrend <- modelIsTrendy & any(parameters=="t");
        parametersCommonSeasonal <- modelIsSeasonal & any(parameters=="s");
        parametersCommonDamped <- modelIsTrendy & damped & any(parameters=="d");
        # Common initials
        initialsCommonLevel <- any(initials=="l");
        initialsCommonTrend <- modelIsTrendy & any(initials=="t");
        initialsCommonSeasonal <- modelIsSeasonal & any(initials=="s");
        # Define common components
        componentsCommonLevel <- any(components=="l");
        componentsCommonTrend <- modelIsTrendy & any(components=="t");
        componentsCommonSeasonal <- modelIsSeasonal & any(components=="s");

        if(componentsCommonLevel){
            componentsCommonTrend[] <- TRUE;
        }

        # Sanity checks. Make initials common if the respective components are
        if(componentsCommonLevel && !initialsCommonLevel){
            initialsCommonLevel[] <- TRUE;
        }
        if(componentsCommonTrend && !initialsCommonTrend){
            initialsCommonTrend[] <- TRUE;
        }
        if(componentsCommonTrend && !parametersCommonDamped){
            parametersCommonDamped[] <- TRUE;
        }
        if(componentsCommonSeasonal && !initialsCommonSeasonal){
            initialsCommonSeasonal[] <- TRUE;
        }

        # Number of parameters in the model
        nParametersLevel <- nSeries^(!parametersCommonLevel);
        nParametersTrend <- modelIsTrendy*nSeries^(!parametersCommonTrend);
        nParametersSeasonal <- modelIsSeasonal*nSeries^(!parametersCommonSeasonal);
        nParametersDamped <- modelIsTrendy*damped*nSeries^(!parametersCommonDamped);

        # Number of overall initials in the model
        nInitialsLevel <- nSeries^(!initialsCommonLevel);
        nInitialsTrend <- modelIsTrendy*nSeries^(!initialsCommonTrend);
        nInitialsSeasonal <- modelIsSeasonal*nSeries^(!initialsCommonSeasonal);

        # Number of overall components in the model
        nComponentsLevel <- nSeries^(!componentsCommonLevel);
        nComponentsTrend <- modelIsTrendy*nSeries^(!componentsCommonTrend);
        nComponentsSeasonal <- modelIsSeasonal*nSeries^(!componentsCommonSeasonal);

        nComponentsNonSeasonal <- nComponentsLevel + nComponentsTrend;
        nComponentsAll <- nComponentsLevel + nComponentsTrend + nComponentsSeasonal;

        # nParamMax <-  nParamMax +
        #     nParametersLevel + nParametersTrend + nParametersSeasonal + nParametersDamped +
        #     nComponentsLevel + nComponentsTrend + nComponentsSeasonal;

        return(list(modelIsTrendy=modelIsTrendy,modelIsSeasonal=modelIsSeasonal, lagsModelMax=lagsModelMax,
                    # PIC booleans
                    parametersCommonLevel=parametersCommonLevel,parametersCommonTrend=parametersCommonTrend,
                    parametersCommonSeasonal=parametersCommonSeasonal,parametersCommonDamped=parametersCommonDamped,
                    initialsCommonLevel=initialsCommonLevel,initialsCommonTrend=initialsCommonTrend,
                    initialsCommonSeasonal=initialsCommonSeasonal,componentsCommonLevel=componentsCommonLevel,
                    componentsCommonTrend=componentsCommonTrend,componentsCommonSeasonal=componentsCommonSeasonal,
                    # PIC number of elements
                    nParametersLevel=nParametersLevel, nParametersTrend=nParametersTrend,
                    nParametersSeasonal=nParametersSeasonal, nParametersDamped=nParametersDamped,
                    nInitialsLevel=nInitialsLevel,nInitialsTrend=nInitialsTrend,
                    nInitialsSeasonal=nInitialsSeasonal,nComponentsLevel=nComponentsLevel,
                    nComponentsTrend=nComponentsTrend,nComponentsSeasonal=nComponentsSeasonal,
                    nComponentsNonSeasonal=nComponentsNonSeasonal,nComponentsAll=nComponentsAll));
    }

    ##### Basic matrices creator #####
    # This thing returns matVt, matF, matG, matW, dampedValue, initialValue
    # and initialSeasonValue if they are not provided + lagsModel
    creatorVETS <- function(...){
        # ParentEnvironment <- list(...)[['ParentEnvironment']];

        #### Blocks for persistence matrix ####
        if(componentsCommonLevel){
            matGBlockAlpha <- matrix(1,1,nSeries);
        }
        else{
            matGBlockAlpha <- diag(nSeries);
        }
        if(modelIsTrendy){
            if(componentsCommonTrend){
                matGBlockBeta <- matrix(1,1,nSeries);
            }
            else{
                matGBlockBeta <- diag(nSeries);
            }
        }
        else{
            matGBlockBeta <- NULL;
        }
        if(modelIsSeasonal){
            if(componentsCommonSeasonal){
                matGBlockGamma <- matrix(1,1,nSeries);
            }
            else{
                matGBlockGamma <- diag(nSeries);
            }
        }
        else{
            matGBlockGamma <- NULL;
        }
        matG <- rbind(matGBlockAlpha, matGBlockBeta, matGBlockGamma);

        if(damped){
            dampedValue <- 0.95;
        }
        else{
            dampedValue <- 1;
        }

        #### Blocks for transition matrix ####
        if(componentsCommonLevel){
            matFBlock11 <- 1;
        }
        else{
            matFBlock11 <- diag(nSeries);
        }
        # L & T elements
        if(modelIsTrendy){
            if(componentsCommonTrend){
                if(componentsCommonLevel){
                    matFBlock12 <- dampedValue;
                    matFBlock21 <- 0;
                }
                else{
                    matFBlock12 <- matrix(dampedValue,nSeries,1);
                    matFBlock21 <- matrix(0,1,nSeries);
                }
                matFBlock22 <- dampedValue;
            }
            else{
                if(componentsCommonLevel){
                    matFBlock12 <- matrix(dampedValue,1,nSeries);
                    matFBlock21 <- matrix(0,nSeries,1);
                }
                else{
                    matFBlock12 <- dampedValue*diag(nSeries);
                    matFBlock21 <- matrix(0,nSeries,nSeries);
                }
                matFBlock22 <- dampedValue*diag(nSeries);
            }
        }
        else{
            matFBlock12 <- NULL;
            matFBlock21 <- NULL;
            matFBlock22 <- NULL;
        }
        # L & S elements
        if(modelIsSeasonal){
            if(componentsCommonSeasonal){
                if(componentsCommonLevel){
                    matFBlock13 <- 0;
                    matFBlock31 <- 0;
                }
                else{
                    matFBlock13 <- matrix(0,nSeries,1);
                    matFBlock31 <- matrix(0,1,nSeries);
                }
                matFBlock33 <- 1;
            }
            else{
                if(componentsCommonLevel){
                    matFBlock13 <- matrix(0,1,nSeries);
                    matFBlock31 <- matrix(0,nSeries,1);
                }
                else{
                    matFBlock13 <- matrix(0,nSeries,nSeries);
                    matFBlock31 <- matrix(0,nSeries,nSeries);
                }
                matFBlock33 <- diag(nSeries);
            }
        }
        else{
            matFBlock13 <- NULL;
            matFBlock31 <- NULL;
            matFBlock33 <- NULL;
        }
        # T & S elements
        if(modelIsTrendy && modelIsSeasonal){
            if(componentsCommonTrend){
                if(componentsCommonSeasonal){
                    matFBlock23 <- 0;
                    matFBlock32 <- 0;
                }
                else{
                    matFBlock23 <- matrix(0,1,nSeries);
                    matFBlock32 <- matrix(0,nSeries,1);
                }
            }
            else{
                if(componentsCommonSeasonal){
                    matFBlock23 <- matrix(0,nSeries,1);
                    matFBlock32 <- matrix(0,1,nSeries);
                }
                else{
                    matFBlock23 <- matrix(0,nSeries,nSeries);
                    matFBlock32 <- matrix(0,nSeries,nSeries);
                }
            }
        }
        else{
            matFBlock23 <- NULL;
            matFBlock32 <- NULL;
        }
        matF <- rbind(cbind(matFBlock11,matFBlock12,matFBlock13,deparse.level=0),
                      cbind(matFBlock21,matFBlock22,matFBlock23,deparse.level=0),
                      cbind(matFBlock31,matFBlock32,matFBlock33,deparse.level=0));

        #### Blocks for measurement matrix ####
        if(componentsCommonLevel){
            matWBlock1 <- matrix(1,nSeries,1);
        }
        else{
            matWBlock1 <- diag(nSeries);
        }
        if(modelIsTrendy){
            if(componentsCommonTrend){
                matWBlock2 <- matrix(dampedValue,nSeries,1);
            }
            else{
                matWBlock2 <- dampedValue*diag(nSeries);
            }
        }
        else{
            matWBlock2 <- NULL;
        }
        if(modelIsSeasonal){
            if(componentsCommonSeasonal){
                matWBlock3 <- matrix(1,nSeries,1);
            }
            else{
                matWBlock3 <- diag(nSeries);
            }
        }
        else{
            matWBlock3 <- NULL;
        }
        matW <- cbind(matWBlock1,matWBlock2,matWBlock3,deparse.level=0);

        #### Vector of states ####
        matVt <- matrix(NA, nComponentsAll, obsStates);

        XValues <- rbind(rep(1,obsInSample),c(1:obsInSample));
        initialValue <- switch(Etype, "M"=log(yInSample), yInSample) %*% t(XValues) %*% solve(XValues %*% t(XValues));

        #### !!! This is a temporary fix for log(0) on intermittent demand ####
        if(any(is.nan(initialValue))){
            initialValue[,1] <- rowMeans(yInSample);
            initialValue[,2] <- switch(Ttype, "M"=1, 0);
        }

        j <- 0;
        # Fill in level initials
        if(initialsCommonLevel){
            matVt[j+1:nComponentsLevel,1:lagsModelMax] <- initialValueNew01 <- mean(initialValue[,1]);
        }
        else{
            matVt[j+1:nComponentsLevel,1:lagsModelMax] <- initialValueNew01 <- initialValue[,1];
        }
        j[] <- j+nComponentsLevel;

        # Fill in trend initials
        if(modelIsTrendy){
            if(initialsCommonTrend){
                matVt[j+1:nComponentsTrend,1:lagsModelMax] <- initialValueNew02 <- mean(initialValue[,2]);
            }
            else{
                matVt[j+1:nComponentsTrend,1:lagsModelMax] <- initialValueNew02 <- initialValue[,2];
            }
            j[] <- j+nComponentsTrend;
        }
        else{
            initialValueNew02 <- NULL;
        }
        initialValue <- c(initialValueNew01, initialValueNew02);

        if(modelIsSeasonal){
            # Matrix of dummies for seasons
            XValues <- matrix(rep(diag(lagsModelMax),ceiling(obsInSample/lagsModelMax)),lagsModelMax)[,1:obsInSample];
            initialSeasonValue <- (switch(Etype, "M"=log(yInSample), yInSample)-
                                       rowMeans(switch(Etype, "M"=log(yInSample), yInSample))) %*%
                t(XValues) %*% solve(XValues %*% t(XValues));
            # Renormalise initials
            initialSeasonValue[] <- initialSeasonValue - rowMeans(initialSeasonValue);

            if(initialsCommonSeasonal){
                matVt[j+1:nComponentsSeasonal,1:lagsModelMax] <- initialSeasonValue <- matrix(colMeans(initialSeasonValue),1);
            }
            else{
                matVt[j+1:nComponentsSeasonal,1:lagsModelMax] <- t(initialSeasonValue);
            }
            # Remove one of columns, to preserve degrees of freedom (normalisation)
            initialSeasonValue <- initialSeasonValue[,-1,drop=FALSE];
        }
        else{
            initialSeasonValue <- NULL;
        }

        #### Names of states ####
        if(componentsCommonLevel){
            statesNames <- "level";
        }
        else{
            statesNames <- paste0("level",c(1:nSeries));
        }
        if(modelIsTrendy){
            if(componentsCommonTrend){
                statesNames <- c(statesNames,"trend");
            }
            else{
                statesNames <- c(statesNames,paste0("trend",c(1:nSeries)));
            }
        }
        if(modelIsSeasonal){
            if(componentsCommonSeasonal){
                statesNames <- c(statesNames,"seasonal");
            }
            else{
                statesNames <- c(statesNames,paste0("seasonal",c(1:nSeries)));
            }
        }
        rownames(matVt) <- statesNames;
        rownames(matF) <- colnames(matF) <- statesNames;
        colnames(matW) <- rownames(matG) <- statesNames;
        rownames(matW) <- colnames(matG) <- paste0("Series",c(1:nSeries));

        ### lagsModel vector
        lagsModel <- matrix(1,nComponentsAll,1);
        if(modelIsSeasonal){
            lagsModel[nComponentsNonSeasonal + 1:nComponentsSeasonal] <- lagsModelMax;
        }

        return(list(matVt=matVt, matF=matF, matG=matG, matW=matW,
                    # initialValues are needed for the initialiser
                    lagsModel=lagsModel, initialValue=initialValue, initialSeasonValue=initialSeasonValue));
    }

    ##### Basic matrices filler #####
    # This thing fills in matVt, matF, matG and matW with values from B and returns the corrected values
    fillerVETS <- function(matVt, matF, matG, matW, B,
                           lagsModelMax, nSeries, modelIsTrendy, modelIsSeasonal, damped,
                           componentsCommonLevel, componentsCommonTrend, componentsCommonSeasonal,
                           nParametersLevel, nParametersTrend, nParametersSeasonal, nParametersDamped,
                           nInitialsLevel, nInitialsTrend, nInitialsSeasonal,
                           nComponentsLevel, nComponentsTrend, nComponentsSeasonal){

        j <- 0;
        k <- 0;
        #### Blocks for persistence matrix ####
        # alpha
        if(componentsCommonLevel){
            matG[1:nComponentsLevel,1:nSeries] <- B[1:nParametersLevel];
        }
        else{
            matG[1:nComponentsLevel,1:nSeries] <- B[1:nParametersLevel]*diag(nSeries);
        }
        j[] <- j+nParametersLevel;
        k[] <- k+nComponentsLevel;
        # beta
        if(modelIsTrendy){
            if(componentsCommonTrend){
                matG[k+1:nComponentsTrend,1:nSeries] <- B[j+1:nParametersTrend];
            }
            else{
                matG[k+1:nComponentsTrend,1:nSeries] <- B[j+1:nParametersTrend]*diag(nSeries);
            }
            j[] <- j+nParametersTrend;
            k[] <- k+nComponentsTrend;
        }
        # gamma
        if(modelIsSeasonal){
            if(componentsCommonSeasonal){
                matG[k+1:nComponentsSeasonal,1:nSeries] <- B[j+1:nParametersSeasonal];
            }
            else{
                matG[k+1:nComponentsSeasonal,1:nSeries] <- B[j+1:nParametersSeasonal]*diag(nSeries);
            }
            j[] <- j+nParametersSeasonal;
        }

        #### Blocks for dampening ####
        if(modelIsTrendy && damped){
            # Transition matrix
            if(componentsCommonTrend){
                if(componentsCommonLevel){
                    matF[1,nComponentsLevel+1] <- B[j+1:nParametersDamped];
                }
                else{
                    matF[1:nComponentsLevel,nComponentsLevel+1] <- B[j+1:nParametersDamped];
                }
                matF[nComponentsLevel+1,nComponentsLevel+1] <- B[j+1:nParametersDamped];
            }
            else{
                if(componentsCommonLevel){
                    matF[1,1+1:nComponentsTrend] <- B[j+1:nParametersDamped];
                }
                else{
                    matF[1:nComponentsLevel,nComponentsLevel+1:nComponentsTrend] <- B[j+1:nParametersDamped]*diag(nComponentsTrend);
                }
                matF[nComponentsLevel+1:nComponentsTrend,nComponentsLevel+1:nComponentsTrend] <-
                    B[j+1:nParametersDamped]*diag(nComponentsTrend);
            }

            # Measurement matrix
            if(componentsCommonTrend){
                matW[1:nSeries,nComponentsLevel+1] <- B[j+1:nParametersDamped];
            }
            else{
                matW[1:nSeries,nComponentsLevel+1:nSeries] <- B[j+1:nParametersDamped]*diag(nComponentsTrend);
            }
            j[] <- j+nParametersDamped;
        }

        k <- 0;
        # Fill in level initials
        matVt[1:nComponentsLevel,1:lagsModelMax] <- B[j+1:nInitialsLevel];
        k[] <- k+nComponentsLevel;
        j[] <- j+nInitialsLevel;

        # Fill in trend initials
        if(modelIsTrendy){
            matVt[k+1:nComponentsTrend,1:lagsModelMax] <- B[j+1:nInitialsTrend];
            k[] <- k+nComponentsTrend;
            j[] <- j+nInitialsTrend;
        }

        if(modelIsSeasonal){
            matVt[k+1:nComponentsSeasonal,2:lagsModelMax] <- matrix(B[j+1:(nInitialsSeasonal*(lagsModelMax-1))],
                                                                    nComponentsSeasonal,(lagsModelMax-1),byrow=TRUE);
            # Normalise initials
            matVt[k+1:nComponentsSeasonal,1] <- -rowSums(matVt[k+1:nComponentsSeasonal,2:lagsModelMax]);

        }

        return(list(matVt=matVt,matF=matF,matG=matG,matW=matW));
    }

    ##### B values for estimation #####
    # Function constructs default bounds where B values should lie
    initialiserVETS <- function(lagsModelMax, nSeries, modelIsTrendy, modelIsSeasonal,
                                componentsCommonLevel, componentsCommonTrend, componentsCommonSeasonal,
                                nComponentsNonSeasonal, nInitialsLevel, nInitialsTrend, nInitialsSeasonal,
                                parametersCommonLevel, parametersCommonTrend, parametersCommonSeasonal, parametersCommonDamped,
                                nParametersLevel, nParametersTrend, nParametersSeasonal, nParametersDamped,
                                initialsCommonLevel, initialsCommonTrend, initialsCommonSeasonal,
                                initialValue, initialSeasonValue){
        # Smoothing parameters, Dampening parameter, initials
        BLower <- BUpper <- B <- rep(NA, nParametersLevel + nParametersTrend + nParametersSeasonal + nParametersDamped +
                                         nInitialsLevel + nInitialsTrend + nInitialsSeasonal*(lagsModelMax-1));

        j <- 0;
        # alpha
        B[1:nParametersLevel] <- 0.1;
        BLower[1:nParametersLevel] <- -5;
        BUpper[1:nParametersLevel] <- 5;
        names(B)[1:nParametersLevel] <- paste0("alpha",c(1:nParametersLevel));
        j[] <- j+nParametersLevel;
        # beta
        if(modelIsTrendy){
            B[j+1:nParametersTrend] <- 0.05;
            BLower[j+1:nParametersTrend] <- -5;
            BUpper[j+1:nParametersTrend] <- 5;
            names(B)[j+1:nParametersTrend] <- paste0("beta",c(1:nParametersTrend));
            j[] <- j+nParametersTrend;
        }
        # gamma
        if(modelIsSeasonal){
            B[j+1:nParametersSeasonal] <- 0.05;
            BLower[j+1:nParametersSeasonal] <- -5;
            BUpper[j+1:nParametersSeasonal] <- 5;
            names(B)[j+1:nParametersSeasonal] <- paste0("gamma",c(1:nParametersSeasonal));
            j[] <- j+nParametersSeasonal;
        }
        # damped
        if(modelIsTrendy && damped){
            B[j+1:nParametersDamped] <- 0.95;
            BLower[j+1:nParametersDamped] <- 0;
            BUpper[j+1:nParametersDamped] <- 1;
            names(B)[j+1:nParametersDamped] <- paste0("phi",c(1:nParametersDamped));
            j[] <- j+nParametersDamped;
        }

        # Initial level and trend
        B[j+1:nInitialsLevel] <- initialValue[1:nInitialsLevel];
        BLower[j+1:nInitialsLevel] <- -Inf;
        BUpper[j+1:nInitialsLevel] <- Inf;
        names(B)[j+1:nInitialsLevel] <- paste0("level",c(1:nInitialsLevel));
        if(modelIsTrendy){
            B[j+nInitialsLevel+1:nInitialsTrend] <- initialValue[nInitialsLevel+1:nInitialsTrend];
            BLower[j+nInitialsLevel+1:nInitialsTrend] <- -Inf;
            BUpper[j+nInitialsLevel+1:nInitialsTrend] <- Inf;
            names(B)[j+nInitialsLevel+1:nInitialsTrend] <- paste0("trend",c(1:nInitialsTrend));
        }
        j[] <- j+nComponentsNonSeasonal;

        # Initial seasonal components
        if(modelIsSeasonal){
            # -1 is due to normalisation of seasonal states
            B[j+1:(nInitialsSeasonal*(lagsModelMax-1))] <- initialSeasonValue;
            BLower[j+1:(nInitialsSeasonal*(lagsModelMax-1))] <- -Inf;
            BUpper[j+1:(nInitialsSeasonal*(lagsModelMax-1))] <- Inf;
            names(B)[j+1:(nInitialsSeasonal*(lagsModelMax-1))] <- paste0(rep(paste0("seasonal",c(1:nInitialsSeasonal),"_"),
                                                                             each=(lagsModelMax-1)),
                                                                         c(1:(lagsModelMax-1)));
            j[] <- j+nInitialsSeasonal*(lagsModelMax-1);
        }

        return(list(B=B,BLower=BLower,BUpper=BUpper));
    }

    ##### Cost Function for VETS #####
    CF <- function(B){
        elements <- fillerVETS(matVt, matF, matG, matW, B,
                               lagsModelMax, nSeries, modelIsTrendy, modelIsSeasonal, damped,
                               componentsCommonLevel, componentsCommonTrend, componentsCommonSeasonal,
                               nParametersLevel, nParametersTrend, nParametersSeasonal, nParametersDamped,
                               nInitialsLevel, nInitialsTrend, nInitialsSeasonal,
                               nComponentsLevel, nComponentsTrend, nComponentsSeasonal);

        # Check the bounds
        if(bounds=="a"){
            eigenValues <- eigen(elements$matF - elements$matG %*% elements$matW, only.values=TRUE, symmetric=TRUE)$values;
            if(max(abs(eigenValues)>(1 + 1E-50))){
                return(max(abs(eigenValues))*1E+100);
            }
        }

        # Fit the model
        fitting <- vFitterWrap(switch(Etype, "M"=log(yInSample), yInSample),
                               elements$matVt, elements$matF, elements$matW, elements$matG,
                               lagsModel, Etype, Ttype, Stype, ot);

        # Calculate the loss
        if(loss=="l"){
            cfRes <- suppressWarnings(log(det((fitting$errors / normalizer) %*% t(fitting$errors / normalizer) / otObs)) +
                                          nSeries * log(normalizer^2));
        }
        else if(loss=="d"){
            cfRes <- sum(log(apply(fitting$errors^2, 2, sum) / obsInSample));
        }
        else{
            cfRes <- sum(apply(fitting$errors^2, 2, sum) / obsInSample);
        }

        # cfRes <- vOptimiserWrap(yInSample, elements$matVt, elements$matF, elements$matW, elements$matG,
        #                         lagsModel, Etype, Ttype, Stype, loss, normalizer, bounds, ot, otObs);
        # multisteps, initialType, bounds,

        if(is.nan(cfRes) | is.na(cfRes) | is.infinite(cfRes)){
            cfRes <- 1e+100;
        }

        return(cfRes);
    }

    ##### Basic estimation function for vets() #####
    estimatorVETS <- function(...){
        environment(creatorVETS) <- environment();
        environment(initialiserVETS) <- environment();
        environment(vLikelihoodFunction) <- environment();
        environment(vICFunction) <- environment();
        environment(CF) <- environment();

        list2env(architectorVETS(Etype, Ttype, Stype, damped, nSeries),environment());
        elements <- creatorVETS();
        list2env(elements,environment());
        modelCurrent <- paste0(Etype, Ttype, ifelse(modelIsTrendy & damped,"d",""), Stype);

        BList <- initialiserVETS(lagsModelMax, nSeries, modelIsTrendy, modelIsSeasonal,
                                 componentsCommonLevel, componentsCommonTrend, componentsCommonSeasonal,
                                 nComponentsNonSeasonal, nInitialsLevel, nInitialsTrend, nInitialsSeasonal,
                                 parametersCommonLevel, parametersCommonTrend, parametersCommonSeasonal, parametersCommonDamped,
                                 nParametersLevel, nParametersTrend, nParametersSeasonal, nParametersDamped,
                                 initialsCommonLevel, initialsCommonTrend, initialsCommonSeasonal,
                                 initialValue, initialSeasonValue);
        if(is.null(B) && is.null(ub) && is.null(lb)){
            B <- BList$B;

            if(any((B>=BList$BUpper),(B<=BList$BLower))){
                B[B>=BList$BUpper] <- BList$BUpper[B>=BList$BUpper] - 0.01;
                B[B<=BList$BLower] <- BList$BLower[B<=BList$BLower] + 0.01;
            }
        }
        else{
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

        # Parameters are chosen to speed up the optimisation process and have decent accuracy
        res <- nloptr(B, CF, lb=BList$BLower, ub=BList$BUpper,
                      opts=list(algorithm=algorithm1, xtol_rel=xtol_rel1, maxeval=maxeval, print_level=print_level));
        B[] <- res$solution;

        # This is just in case something went out of the bounds
        if(any((B>=BList$BUpper),(B<=BList$BLower))){
            BList$BUpper[B>=BList$BUpper] <- B[B>=BList$BUpper] + 1;
            BList$BLower[B<=BList$BLower] <- B[B<=BList$BLower] - 1;
        }

        if(print_level_hidden>0){
            print(res);
        }

        res2 <- nloptr(B, CF, lb=BList$BLower, ub=BList$BUpper,
                       opts=list(algorithm=algorithm2, xtol_rel=xtol_rel2, maxeval=maxeval, print_level=print_level));
        # This condition is needed in order to make sure that we did not make the solution worse
        if(res2$objective <= res$objective){
            res <- res2;
        }
        B[] <- res$solution;

        if(print_level_hidden>0){
            print(res);
        }

        if(all(B==BList$B)){
            warning(paste0("Failed to optimise the model ETS(", modelCurrent,
                           "). Try different initialisation maybe?\nAnd check all the messages and warnings...",
                           "If you did your best, but the optimiser still fails, report this to the maintainer, please."),
                    call.=FALSE);
        }

        # First part is for the covariance matrix
        if(loss=="l"){
            nParam <- nSeries * (nSeries + 1) / 2 + length(B);
        }
        else{
            nParam <- nSeries + length(B);
        }

        ICValues <- vICFunction(nParam=nParam,B=B,Etype=Etype);
        ICs <- ICValues$ICs;
        logLik <- ICValues$llikelihood;

        # Write down Fisher Information if needed
        if(FI){
            environment(vLikelihoodFunction) <- environment();
            FI <- -numDeriv::hessian(vLikelihoodFunction,B);
            rownames(FI) <- BList$BNames;
            colnames(FI) <- BList$BNames;
        }

        return(list(ICs=ICs,objective=res$objective,B=B,nParam=nParam,logLik=logLik,FI=FI));
    }

    ##### Function selects ETS components #####
    selectorVETS <- function(silent=TRUE,...){
        environment(estimatorVETS) <- environment();

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
        vetsModels <- vector("list", length(modelsPool));
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
            vetsModels[[i]] <- estimatorVETS(ParentEnvironment=environment());
        }
        # Prepare the return of the best model
        vetsModelsICs <- sapply(vetsModels,"[[","ICs");
        colnames(vetsModelsICs) <- modelsPool;
        iBest <- which.min(vetsModelsICs[ic,]);
        vetsModels[[iBest]]$model <- modelsPool[iBest];
        vetsModels[[iBest]]$Etype <- substring(modelsPool[iBest],1,1);
        vetsModels[[iBest]]$Ttype <- substring(modelsPool[iBest],2,2);
        if(nchar(modelsPool[iBest])==4){
            vetsModels[[iBest]]$damped <- TRUE;
            vetsModels[[iBest]]$Stype <- substring(modelsPool[iBest],4,4);
        }
        else{
            vetsModels[[iBest]]$damped <- FALSE;
            vetsModels[[iBest]]$Stype <- substring(modelsPool[iBest],3,3);
        }
        vetsModels[[iBest]]$ICsAll <- vetsModelsICs;
        vetsModels[[iBest]]$ICs <- vetsModelsICs[,iBest];
        vetsModels[[iBest]]$icBest <- vetsModelsICs[ic,iBest];
        # Rename "objective" into "cfObjective"
        names(vetsModels[[iBest]])[names(vetsModels[[iBest]])=="objective"] <- "cfObjective";
        return(vetsModels[[iBest]]);
    }

    ##### Function constructs the VETS function #####
    callerVETS <- function(silent=FALSE,...){
        if(modelDo=="estimate"){
            environment(estimatorVETS) <- environment();
            res <- estimatorVETS(ParentEnvironment=environment());
            listToReturn <- list(Etype=Etype,Ttype=Ttype,Stype=Stype,damped=damped,
                                 cfObjective=res$objective,B=res$B,ICs=res$ICs,icBest=res$ICs[ic],
                                 ICsAll=res$ICs,nParam=res$nParam,logLik=res$logLik,FI=res$FI);

            return(listToReturn);
        }
        else if(modelDo=="select"){
            return(selectorVETS(ParentEnvironment=environment()));
        }
        # else{
        #     environment(CF) <- environment();
        #     environment(vICFunction) <- environment();
        #     environment(vLikelihoodFunction) <- environment();
        #     environment(creatorVETS) <- environment();
        #     elements <- creatorVETS();
        #     list2env(elements,environment());
        #
        #     B <- c(persistenceValue);
        #     BNames <- paste0("Persistence",c(1:length(persistenceValue)));
        #     if(damped){
        #         B <- c(B,dampedValue);
        #         BNames <- c(BNames,paste0("phi",c(1:length(dampedValue))));
        #     }
        #     if(transitionType=="d"){
        #         transitionLength <- length(B);
        #         # Write values from the rest of transition matrix
        #         for(i in 1:nSeries){
        #             B <- c(B, c(transitionValue[c(1:nComponentsAll)+nComponentsAll*(i-1),
        #                                         setdiff(c(1:nSeries*nComponentsAll),
        #                                                 c(1:nComponentsAll)+nComponentsAll*(i-1))]));
        #         }
        #         transitionLength <- length(B) - transitionLength;
        #         BNames <- c(BNames,paste0("transition",c(1:transitionLength)));
        #     }
        #     B <- c(B,initialValue);
        #     BNames <- c(BNames,paste0("initial",c(1:length(initialValue))));
        #     if(Stype!="N"){
        #         B <- c(B,initialSeasonValue);
        #         BNames <- c(BNames,paste0("initialSeason",c(1:length(initialSeasonValue))));
        #     }
        #     names(B) <- BNames;
        #
        #     cfObjective <- CF(B);
        #
        #     # Number of parameters
        #     # First part is for the covariance matrix
        #     if(loss=="l"){
        #         nParam <- nSeries * (nSeries + 1) / 2;
        #     }
        #     else if(loss=="d"){
        #         nParam <- nSeries;
        #     }
        #     else{
        #         nParam <- nSeries;
        #     }
        #
        #     ICValues <- vICFunction(nParam=nParam,B=B,Etype=Etype);
        #     logLik <- ICValues$llikelihood;
        #     ICs <- ICValues$ICs;
        #     icBest <- ICs[ic];
        #
        #     # Write down Fisher Information if needed
        #     if(FI){
        #         environment(vLikelihoodFunction) <- environment();
        #         FI <- -numDeriv::hessian(vLikelihoodFunction,B);
        #         rownames(FI) <- BNames;
        #         colnames(FI) <- BNames;
        #     }
        #
        #     listToReturn <- list(Etype=Etype,Ttype=Ttype,Stype=Stype,damped=damped,
        #                          cfObjective=cfObjective,B=B,ICs=ICs,icBest=icBest,
        #                          nParam=nParam,logLik=logLik,FI=FI);
        #     return(listToReturn);
        # }
    }

    ##### Preset yFitted, yForecast, errors and basic parameters #####
    yFitted <- matrix(NA,nSeries,obsInSample);
    yForecast <- matrix(NA,nSeries,h);
    errors <- matrix(NA,nSeries,obsInSample);
    rownames(yFitted) <- rownames(yForecast) <- rownames(errors) <- dataNames;

    ##### Now do estimation and model selection #####
    environment(creatorVETS) <- environment();
    environment(fillerVETS) <- environment();
    environment(vssFitter) <- environment();
    environment(vssForecaster) <- environment();

    ##### Fit the model and produce forecast #####
    list2env(callerVETS(silent=silentText),environment());
    list2env(architectorVETS(Etype, Ttype, Stype, damped, nSeries),environment());
    list2env(creatorVETS(),environment());
    list2env(fillerVETS(matVt, matF, matG, matW, B,
                        lagsModelMax, nSeries, modelIsTrendy, modelIsSeasonal, damped,
                        componentsCommonLevel, componentsCommonTrend, componentsCommonSeasonal,
                        nParametersLevel, nParametersTrend, nParametersSeasonal, nParametersDamped,
                        nInitialsLevel, nInitialsTrend, nInitialsSeasonal,
                        nComponentsLevel, nComponentsTrend, nComponentsSeasonal),environment());

    if(Etype=="M"){
        cfObjective <- exp(cfObjective);
    }

    if(damped){
        model <- paste0(Etype,Ttype,"d",Stype);
    }
    else{
        model <- paste0(Etype,Ttype,Stype);
    }

    vssFitter(ParentEnvironment=environment());
    vssForecaster(ParentEnvironment=environment());

    parametersNumber[1,1] <- length(B);

    # This is needed anyway for the reusability of the model
    # if(seasonalType=="i"){
    #     initialPlaces <- nComponentsAll*(c(1:nSeries)-1)+1;
    #     initialNames <- "level";
    #     if(Ttype!="N"){
    #         initialPlaces <- c(initialPlaces,nComponentsAll*(c(1:nSeries)-1)+2);
    #         initialPlaces <- sort(initialPlaces);
    #         initialNames <- c(initialNames,"trend");
    #     }
    #     if(initialEstimate){
    #         initialValue <- matrix(matVt[initialPlaces,lagsModelMax],nComponentsNonSeasonal*nSeries,1);
    #         parametersNumber[1,1] <- parametersNumber[1,1] + length(unique(as.vector(initialValue)));
    #     }
    # }
    # else{
    #     initialNames <- "level";
    #     if(Ttype!="N"){
    #         initialNames <- c(initialNames,"trend");
    #     }
    #     if(initialEstimate){
    #         initialValue <- matrix(matVt[1:(nComponentsNonSeasonal*nSeries),lagsModelMax],
    #                                nComponentsNonSeasonal*nSeries,1);
    #         parametersNumber[1,1] <- parametersNumber[1,1] + length(unique(as.vector(initialValue)));
    #     }
    # }
    # rownames(initialValue) <- paste0(rep(dataNames,each=nComponentsNonSeasonal), "_", initialNames);
    #
    # if(modelIsSeasonal){
    #     if(seasonalType=="i"){
    #         if(initialSeasonEstimate){
    #             initialPlaces <- nComponentsAll*(c(1:nSeries)-1)+nComponentsAll;
    #             initialSeasonValue <- matrix(matVt[initialPlaces,1:lagsModelMax],nSeries,lagsModelMax);
    #             parametersNumber[1,1] <- parametersNumber[1,1] + length(unique(as.vector(initialSeasonValue)));
    #         }
    #         rownames(initialSeasonValue) <- dataNames;
    #     }
    #     else{
    #         initialSeasonValue <- matrix(matVt[nComponentsNonSeasonal*nSeries+1,1:lagsModelMax],1,lagsModelMax);
    #         parametersNumber[1,1] <- parametersNumber[1,1] + lagsModelMax;
    #         rownames(initialSeasonValue) <- "Common";
    #     }
    #     colnames(initialSeasonValue) <- paste0("Seasonal",c(1:lagsModelMax));
    # }

    yFitted <- ts(t(yFitted),start=dataStart,frequency=dataFreq);
    errors <- ts(t(errors),start=dataStart,frequency=dataFreq);

    yForecast <- ts(t(yForecast),start=yForecastStart,frequency=dataFreq);
    if(!is.matrix(yForecast)){
        yForecast <- as.matrix(yForecast,h,nSeries);
    }
    colnames(yForecast) <- dataNames;
    yForecastStart <- start(yForecast)
    if(any(intervalType==c("i","u","l"))){
        PI <-  ts(PI,start=yForecastStart,frequency=dataFreq);
    }

    if(loss=="l"){
        loss <- "likelihood";
        parametersNumber[1,1] <- parametersNumber[1,1] + nSeries * (nSeries + 1) / 2;
    }
    else if(loss=="d"){
        loss <- "diagonal";
        parametersNumber[1,1] <- parametersNumber[1,1] + nSeries;
    }
    else{
        loss <- "trace";
        parametersNumber[1,1] <- parametersNumber[1,1] + nSeries;
    }

    parametersNumber[1,4] <- sum(parametersNumber[1,1:3]);
    parametersNumber[2,4] <- sum(parametersNumber[2,1:3]);

    ##### Now let's deal with the holdout #####
    if(holdout){
        # if(modelIsMultiplicative){
        #     yInSample[] <- exp(yInSample);
        # }
        yHoldout <- ts(y[(obsInSample+1):obsAll,],start=yForecastStart,frequency=dataFreq);
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
        yHoldout <- NA;
        errorMeasures <- NA;
    }

    modelName <- "VETS";
    modelName <- paste0(modelName,"(",model,")");

    if(occurrence!="n"){
        modelName <- paste0("i",modelName);
    }

    # If something is restricted, add "PIC
    if(!all(parametersCommonLevel,parametersCommonTrend,parametersCommonSeasonal,damped & parametersCommonDamped,
            initialsCommonLevel,initialsCommonTrend,initialsCommonSeasonal,
            componentsCommonLevel,componentsCommonTrend,componentsCommonSeasonal)){
        modelName <- paste0(modelName,"PIC");
        submodelNameP <- paste0(ifelse(c(parametersCommonLevel,parametersCommonTrend,
                                         parametersCommonSeasonal,damped & parametersCommonDamped),
                                       c("L","T","S","D"),""),collapse="");
        submodelNameI <- paste0(ifelse(c(initialsCommonLevel,initialsCommonTrend,
                                         initialsCommonSeasonal),
                                       c("L","T","S"),""),collapse="");
        submodelNameC <- paste0(ifelse(c(componentsCommonLevel,componentsCommonTrend,
                                         componentsCommonSeasonal),
                                       c("L","T","S"),""),collapse="");
        if(nchar(submodelNameP)==0){
            submodelNameP <- "N";
        }
        if(nchar(submodelNameI)==0){
            submodelNameI <- "N";
        }
        if(nchar(submodelNameC)==0){
            submodelNameC <- "N";
        }
        modelName <- paste0(modelName,"(",paste(submodelNameP,submodelNameI,submodelNameC,sep=","),")");
    }

    #     if(modelIsSeasonal){
    #         submodelName <- "[";
    #         if(seasonalType=="c"){
    #             submodelName[] <- paste0(submodelName,"C");
    #         }
    #         else{
    #             submodelName[] <- paste0(submodelName,"I");
    #         }
    #
    #         if(persistenceType=="i"){
    #             submodelName[] <- paste0(submodelName,"I");
    #         }
    #         else if(persistenceType=="c"){
    #             submodelName[] <- paste0(submodelName,"CA");
    #         }
    #         else if(persistenceType=="s"){
    #             submodelName[] <- paste0(submodelName,"CS");
    #         }
    #         else if(persistenceType=="d"){
    #             submodelName[] <- paste0(submodelName,"D");
    #         }
    #
    #         if(initialSeasonType=="i"){
    #             submodelName[] <- paste0(submodelName,"I");
    #         }
    #         else{
    #             submodelName[] <- paste0(submodelName,"C");
    #         }
    #         submodelName[] <- paste0(submodelName,"]");
    #         modelName[] <- paste0(modelName,submodelName);
    #     }

    ##### Print output #####
    if(!silentText){
        if(any(abs(eigen(matF - matG %*% matW)$values)>(1 + 1E-10))){
            warning(paste0("Model VETS(",model,") is unstable! ",
                           "Use a different value of 'bounds' parameter to address this issue!"),
                    call.=FALSE);
        }
    }

    ##### Make a plot #####
    # This is a temporary solution
    if(!silentGraph){
        pages <- ceiling(nSeries / 5);
        perPage <- ceiling(nSeries / pages);
        packs <- c(seq(1, nSeries+1, perPage));
        if(packs[length(packs)]<nSeries+1){
            packs <- c(packs,nSeries+1);
        }
        parDefault <- par(no.readonly=TRUE);
        for(j in 1:pages){
            par(mar=c(4,4,2,1),mfcol=c(perPage,1));
            for(i in packs[j]:(packs[j+1]-1)){
                if(any(intervalType==c("u","i","l"))){
                    plotRange <- range(min(y[,i],yForecast[,i],yFitted[,i],PI[,i*2-1]),
                                       max(y[,i],yForecast[,i],yFitted[,i],PI[,i*2]));
                }
                else{
                    plotRange <- range(min(y[,i],yForecast[,i],yFitted[,i]),
                                       max(y[,i],yForecast[,i],yFitted[,i]));
                }
                plot(y[,i],main=paste0(modelName," on ",dataNames[i]),ylab="Y",
                     ylim=plotRange, xlim=range(time(y[,i])[1],time(yForecast)[max(h,1)]),
                     type="l");
                lines(yFitted[,i],col="purple",lwd=2,lty=2);
                if(h>1){
                    if(any(intervalType==c("u","i","l"))){
                        lines(PI[,i*2-1],col="darkgrey",lwd=3,lty=2);
                        lines(PI[,i*2],col="darkgrey",lwd=3,lty=2);

                        polygon(c(seq(dataDeltat*(yForecastStart[2]-1)+yForecastStart[1],
                                      dataDeltat*(end(yForecast)[2]-1)+end(yForecast)[1],dataDeltat),
                                  rev(seq(dataDeltat*(yForecastStart[2]-1)+yForecastStart[1],
                                          dataDeltat*(end(yForecast)[2]-1)+end(yForecast)[1],dataDeltat))),
                                c(as.vector(PI[,i*2]), rev(as.vector(PI[,i*2-1]))), col="lightgray",
                                border=NA, density=10);
                    }
                    lines(yForecast[,i],col="blue",lwd=2);
                }
                else{
                    if(any(intervalType==c("u","i","l"))){
                        points(PI[,i*2-1],col="darkgrey",lwd=3,pch=4);
                        points(PI[,i*2],col="darkgrey",lwd=3,pch=4);
                    }
                    points(yForecast[,i],col="blue",lwd=2,pch=4);
                }
                abline(v=dataDeltat*(yForecastStart[2]-2)+yForecastStart[1],col="red",lwd=2);
            }
        }
        par(parDefault);
    }

    ##### Return values #####
    model <- list(model=modelName,timeElapsed=Sys.time()-startTime,
                  states=ts(t(matVt),start=(time(y)[1] - dataDeltat*lagsModelMax),frequency=dataFreq),
                  persistence=matG, transition=matF, measurement=matW, B=B,
                  # initialType=initialType,initial=initialValue,initialSeason=initialSeasonValue,
                  nParam=parametersNumber, occurrence=ovesModel,
                  y=y,fitted=yFitted,holdout=yHoldout,residuals=errors,Sigma=Sigma,
                  forecast=yForecast,PI=PI,interval=intervalType,level=level,
                  ICs=ICs,ICsAll=ICsAll,logLik=logLik,lossValue=cfObjective,loss=loss,accuracy=errorMeasures,
                  FI=FI);
    return(structure(model,class=c("legion","smooth")));
}
