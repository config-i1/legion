utils::globalVariables(c("nParamMax","nComponentsAll","nComponentsNonSeasonal","nSeries","modelIsSeasonal",
                         "obsInSample","obsAll","lagsModel","persistenceEstimate","persistenceType",
                         "persistenceValue","damped","dampedEstimate","dampedType","transitionType",
                         "initialEstimate","initialSeasonEstimate","initialSeasonValue","initialSeasonType",
                         "modelIsMultiplicative","matG","matW","B","ub","lb", "maxeval", "algorithm1",
                         "algorithm2", "xtol_rel1", "xtol_rel2", "Sigma","yFitted","PI","dataDeltat",
                         "dataFreq","dataStart","otObs","dataNames","seasonalType",
                         "CF","Etype","FI","ICs","Stype","Ttype","cumulative","errors","h","holdout",
                         "initial","initialType","interval","intervalType","is.vsmooth.sim","lagsModelMax",
                         "level","matF","matVt","measures","nParam","normalizer","obsStates","ot",
                         "silentGraph","silentText","transition","transitionEstimate","yInSample"));

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
#' \code{ZZZ} means that the model will be selected based on the chosen
#' information criteria type.
#' ATTENTION! ONLY PURE ADDITIVE AND PURE MULTIPLICATIVE MODELS ARE
#' AVAILABLE!
#' Pure multiplicative models are done as additive model applied to log(data).
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
# vets(Y,model="ANN",h=10,holdout=TRUE)
#'
#' # Multiplicative damped trend model with common parameters
#' # and initial seasonal indices
# vets(Y,model="MMdM",h=10,holdout=TRUE,parameters=c("l","t","s","d"),
#      initials="seasonal")
#'
#'
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

##### Basic VETS architector
### This function will accept Etype, Ttype, Stype and damped and would return:
# nComponentsNonSeasonal, nComponentsAll, lagsModelMax, modelIsSeasonal, modelIsTrendy, obsStates
# This is needed for model selection
# architectorVETS <- function(...){}

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
    initialValue <- yInSample %*% t(XValues) %*% solve(XValues %*% t(XValues));
    j <- 0;
    # Fill in level initials
    if(componentsCommonLevel){
        matVt[j+1:nComponentsLevel,1:lagsModelMax] <- intiialValueNew01 <- mean(initialValue[,1]);
    }
    else{
        matVt[j+1:nComponentsLevel,1:lagsModelMax] <- intiialValueNew01 <- initialValue[,1];
    }
    j[] <- j+nComponentsLevel;

    # Fill in trend initials
    if(modelIsTrendy){
        if(componentsCommonTrend){
            matVt[j+1:nComponentsTrend,1:lagsModelMax] <- intiialValueNew02 <- mean(initialValue[,2]);
        }
        else{
            matVt[j+1:nComponentsTrend,1:lagsModelMax] <- intiialValueNew02 <- initialValue[,2];
        }
        j[] <- j+nComponentsTrend;
    }
    else{
        intiialValueNew02 <- NULL;
    }
    initialValue <- c(intiialValueNew01, intiialValueNew02);

    if(modelIsSeasonal){
        # Matrix of dummies for seasons
        XValues <- matrix(rep(diag(lagsModelMax),ceiling(obsInSample/lagsModelMax)),lagsModelMax)[,1:obsInSample];
        initialSeasonValue <- (yInSample-rowMeans(yInSample)) %*% t(XValues) %*% solve(XValues %*% t(XValues));

        if(componentsCommonSeasonal){
            matVt[j+1:nComponentsSeasonal,1:lagsModelMax] <- initialSeasonValue <- colMeans(initialSeasonValue);
        }
        else{
            matVt[j+1:nComponentsSeasonal,1:lagsModelMax] <- t(initialSeasonValue);
        }
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
    if(componentsCommonSeasonal){
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
    matVt[1:nComponentsLevel,1:lagsModelMax] <- B[j+1:nComponentsLevel];
    k[] <- k+nComponentsLevel;
    j[] <- j+nComponentsLevel;

    # Fill in trend initials
    if(modelIsTrendy){
        matVt[k+1:nComponentsTrend,1:lagsModelMax] <- B[j+1:nComponentsTrend];
        k[] <- k+nComponentsTrend;
        j[] <- j+nComponentsTrend;
    }

    if(modelIsSeasonal){
        matVt[k+1:nComponentsSeasonal,1:lagsModelMax] <- matrix(B[j+1:(nComponentsSeasonal*lagsModelMax)],
                                                                nComponentsSeasonal,lagsModelMax,byrow=TRUE);
    }

    return(list(matVt=matVt,matF=matF,matG=matG,matW=matW));
}

##### B values for estimation #####
# Function constructs default bounds where B values should lie
initialiserVETS <- function(lagsModelMax, nSeries, modelIsTrendy, modelIsSeasonal,
                            componentsCommonLevel, componentsCommonTrend, componentsCommonSeasonal,
                            nComponentsAll, nComponentsNonSeasonal, nComponentsLevel, nComponentsTrend, nComponentsSeasonal,
                            parametersCommonLevel, parametersCommonTrend, parametersCommonSeasonal, parametersCommonDamped,
                            nParametersLevel, nParametersTrend, nParametersSeasonal, nParametersDamped,
                            initialsCommonLevel, initialsCommonTrend, initialsCommonSeasonal,
                            initialValue, initialSeasonValue){
    # Smoothing parameters, Dampening parameter, initials
    BLower <- BUpper <- B <- rep(NA, nParametersLevel + nParametersTrend + nParametersSeasonal + nParametersDamped +
                                     nComponentsNonSeasonal + nComponentsSeasonal*lagsModelMax);

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
    B[j+1:nComponentsNonSeasonal] <- initialValue;
    BLower[j+1:nComponentsNonSeasonal] <- -Inf;
    BUpper[j+1:nComponentsNonSeasonal] <- Inf;
    names(B)[j+1:nComponentsLevel] <- paste0("level",c(1:nComponentsLevel));
    if(modelIsTrendy){
        names(B)[j+nComponentsLevel+1:nComponentsTrend] <- paste0("trend",c(1:nComponentsTrend));
    }
    j[] <- j+nComponentsNonSeasonal;

    # Initial seasonal components
    if(modelIsSeasonal){
        B[j+1:(nComponentsSeasonal*lagsModelMax)] <- initialSeasonValue;
        BLower[j+1:(nComponentsSeasonal*lagsModelMax)] <- -Inf;
        BUpper[j+1:(nComponentsSeasonal*lagsModelMax)] <- Inf;
        names(B)[j+1:(nComponentsSeasonal*lagsModelMax)] <- paste0(rep(paste0("seasonal",c(1:nComponentsSeasonal),"_"),
                                                                       each=lagsModelMax),
                                                                   c(1:lagsModelMax));
        j[] <- j+nComponentsSeasonal*lagsModelMax;
    }

    return(list(B=B,BLower=BLower,BUpper=BUpper));
}

##### Cost Function for VETS #####
CF <- function(B){
    elements <- fillerVETS(matVt, matF, matG, matW, B,
                           lagsModelMax, nSeries, modelIsTrendy, modelIsSeasonal, damped,
                           componentsCommonLevel, componentsCommonTrend, componentsCommonSeasonal,
                           nParametersLevel, nParametersTrend, nParametersSeasonal, nParametersDamped,
                           nComponentsLevel, nComponentsTrend, nComponentsSeasonal);

    # Check the bounds
    if(bounds=="a"){
        eigenValues <- eigen(elements$matF - elements$matG %*% elements$matW, only.values=TRUE, symmetric=TRUE)$values;
        if(max(abs(eigenValues)>(1 + 1E-50))){
            return(max(abs(eigenValues))*1E+100);
        }
    }

    # Fit the model
    fitting <- vFitterWrap(yInSample, elements$matVt, elements$matF, elements$matW, elements$matG,
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
EstimatorVETS <- function(...){
    environment(creatorVETS) <- environment();
    environment(initialiserVETS) <- environment();
    environment(vLikelihoodFunction) <- environment();
    environment(vICFunction) <- environment();
    environment(CF) <- environment();
    elements <- creatorVETS();
    list2env(elements,environment());

    BList <- initialiserVETS(lagsModelMax, nSeries, modelIsTrendy, modelIsSeasonal,
                             componentsCommonLevel, componentsCommonTrend, componentsCommonSeasonal,
                             nComponentsAll, nComponentsNonSeasonal, nComponentsLevel, nComponentsTrend, nComponentsSeasonal,
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

    normalizer <- sum(colMeans(abs(diff(t(yInSample))),na.rm=TRUE));

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
    B <- res$solution;

    if(print_level_hidden>0){
        print(res);
    }

    if(all(B==BList$B) & modelDo=="estimate"){
        warning(paste0("Failed to optimise the model ETS(", modelCurrent,
                       "). Try different initialisation maybe?\nAnd check all the messages and warnings...",
                       "If you did your best, but the optimiser still fails, report this to the maintainer, please."),
                call.=FALSE);
    }
    names(B) <- BList$BNames;

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

##### Function constructs the VETS function #####
CreatorVETS <- function(silent=FALSE,...){
    if(modelDo=="estimate"){
        environment(EstimatorVETS) <- environment();
        res <- EstimatorVETS(ParentEnvironment=environment());
        listToReturn <- list(Etype=Etype,Ttype=Ttype,Stype=Stype,damped=damped,
                             cfObjective=res$objective,B=res$B,ICs=res$ICs,icBest=res$ICs[ic],
                             nParam=res$nParam,logLik=res$logLik,FI=res$FI);

        return(listToReturn);
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

##### Define modelDo #####
    # if(any(persistenceEstimate, transitionEstimate, dampedEstimate, initialEstimate, initialSeasonEstimate)){
        modelDo <- "estimate";
        modelCurrent <- model;
    # }
    # else{
    #     modelDo <- "nothing";
    #     modelCurrent <- model;
    #     bounds <- "n";
    # }

##### Now do estimation and model selection #####
    environment(creatorVETS) <- environment();
    environment(fillerVETS) <- environment();
    environment(vssFitter) <- environment();
    environment(vssForecaster) <- environment();

    vetsValues <- CreatorVETS(silent=silentText);

##### Fit the model and produce forecast #####
    list2env(vetsValues,environment());
    list2env(creatorVETS(),environment());
    list2env(fillerVETS(matVt, matF, matG, matW, B,
                        lagsModelMax, nSeries, modelIsTrendy, modelIsSeasonal, damped,
                        componentsCommonLevel, componentsCommonTrend, componentsCommonSeasonal,
                        nParametersLevel, nParametersTrend, nParametersSeasonal, nParametersDamped,
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
        if(modelIsMultiplicative){
            yInSample[] <- exp(yInSample);
        }
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
    if(!all(parametersCommonLevel,parametersCommonTrend,parametersCommonSeasonal,parametersCommonDamped,
            initialsCommonLevel,initialsCommonTrend,initialsCommonSeasonal,
            componentsCommonLevel,componentsCommonTrend,componentsCommonSeasonal)){
        modelName <- paste0(modelName,"PIC");
        submodelName <- paste(paste0(ifelse(c(parametersCommonLevel,parametersCommonTrend,
                                               parametersCommonSeasonal,parametersCommonDamped),
                                             c("L","T","S","D"),""),collapse=""),
                               paste0(ifelse(c(initialsCommonLevel,initialsCommonTrend,
                                               initialsCommonSeasonal),
                                             c("L","T","S"),""),collapse=""),
                               paste0(ifelse(c(componentsCommonLevel,componentsCommonTrend,
                                               componentsCommonSeasonal),
                                             c("L","T","S"),""),collapse=""),sep=",");
        modelName <- paste0(modelName,"(",submodelName,")");
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
                  ICs=ICs,logLik=logLik,lossValue=cfObjective,loss=loss,accuracy=errorMeasures,
                  FI=FI);
    return(structure(model,class=c("legion","smooth")));
}
