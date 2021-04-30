utils::globalVariables(c("initialSeason","persistence","phi","otObs",
                         "occurrence","ovesModel","occurrenceModelProvided","seasonal","lags",
                         "loss"));

##### *Checker of input of vector functions* #####
vssInput <- function(smoothType=c("ves","vets"),ParentEnvironment,...){
    smoothType <- smoothType[1];

    ellipsis <- list(...);

    ##### silent #####
    # Fix for cases with TRUE/FALSE.
    if(!is.logical(silent)){
        warning("The parameter silent can only be TRUE or FALSE. Switching to TRUE.",
                call.=FALSE);
        silent <- TRUE;
    }

    #### Check horizon ####
    if(h<=0){
        warning(paste0("You have set forecast horizon equal to ",h,". We hope you know, what you ",
                       "are doing."),
                call.=FALSE);
        if(h<0){
            warning("And by the way, we can't do anything with negative horizon, so we will set ",
                    "it equal to zero.",
                    call.=FALSE);
            h <- 0;
        }
    }

    #### Check data ####
    if(any(is.legion.sim(data))){
        data <- data$data;
        if(length(dim(data))==3){
            warning("Simulated data contains several samples. Selecting a random one.",call.=FALSE);
            data <- ts(data[,,runif(1,1,dim(data)[3])]);
        }
    }

    if(!is.data.frame(data) && !is.numeric(data)){
        stop("The provided data is not a numeric matrix! Can't construct any model!", call.=FALSE);
    }

    if(is.null(dim(data))){
        stop("The provided data is not a matrix or a data.frame! If it is a vector, please use ",
             "es() function instead.",
             call.=FALSE);
    }

    # Get indices and classes of data
    yIndex <- try(time(data),silent=TRUE);
    # If we cannot extract time, do something
    if(inherits(yIndex,"try-error")){
        if(!is.data.frame(data) && !is.null(dim(data))){
            yIndex <- as.POSIXct(rownames(data));
        }
        else if(is.data.frame(data)){
            yIndex <- c(1:nrow(data));
        }
        else{
            yIndex <- c(1:length(data));
        }
    }
    yClasses <- class(data);

    # If this is just a numeric variable, use ts class
    if(all(yClasses=="integer") || all(yClasses=="numeric") ||
       all(yClasses=="data.frame") || all(yClasses=="matrix")){
        if(any(class(yIndex) %in% c("POSIXct","Date"))){
            yClasses <- "zoo";
        }
        else{
            yClasses <- "ts";
        }
    }

    if(is.data.frame(data)){
        data <- as.matrix(data);
    }

    # Number of series in the matrix
    nSeries <- ncol(data);

    correlatedSeries <- cor(data)[upper.tri(cor(data))];
    if(any(correlatedSeries>0.999)){
        warning("Some of series are almost perfectly correlated. This might cause difficulties ",
                "in the estimation. ",
                "Please, try removing some of them if you encounter any problems.",
                call.=FALSE);
    }

    if(is.null(ncol(data))){
        stop("The provided data is not a matrix! Use es() or adam() function instead!",
             call.=FALSE);
    }
    if(ncol(data)==1){
        stop("The provided data contains only one column. Use es() or adam() function instead!",
             call.=FALSE);
    }
    # Check the data for NAs
    if(any(is.na(data))){
        if(!silent){
            warning("Data contains NAs. These observations will be substituted by zeroes.",
                    call.=FALSE);
        }
        data[is.na(data)] <- 0;
    }

    # Define obs, the number of observations of in-sample
    obsInSample <- nrow(data) - holdout*h;

    # Define obsAll, the overal number of observations (in-sample + holdout)
    obsAll <- nrow(data) + (1 - holdout)*h;

    # If obsInSample is negative, this means that we can't do anything...
    if(obsInSample<=0){
        stop("Not enough observations in sample.",
             call.=FALSE);
    }
    # Define the actual values. Transpose the matrix!
    yInSample <- t(data[1:obsInSample,,drop=FALSE]);
    #### For now we just get the first lag. This would need to be modified for multiple seasonal models
    lags <- unique(lags);
    lagsModel <- lags;
    yFrequency <- frequency(data);
    yStart <- yIndex[1];
    yDeltat <- yIndex[2]-yIndex[1];
    if(holdout){
        yForecastStart <- yIndex[obsInSample+1];
        yForecastIndex <- yIndex[-c(1:obsInSample)];
        yInSampleIndex <- yIndex[c(1:obsInSample)];
        yIndexAll <- yIndex;
    }
    else{
        yInSampleIndex <- yIndex;
        yIndexDiff <- diff(tail(yIndex,2));
        yForecastStart <- yIndex[obsInSample]+yIndexDiff;
        if(any(yClasses=="ts")){
            yForecastIndex <- yIndex[obsInSample]+as.numeric(yIndexDiff)*c(1:max(h,1));
        }
        else{
            yForecastIndex <- yIndex[obsInSample]+yIndexDiff*c(1:max(h,1));
        }
        yIndexAll <- c(yIndex,yForecastIndex);
    }

    dataNames <- colnames(data);
    if(!is.null(dataNames)){
        dataNames <- make.names(dataNames);
    }
    else{
        dataNames <- paste0("Series",c(1:nSeries));
    }

    # Number of parameters to estimate / provided
    parametersNumber <- matrix(0,2,4,
                               dimnames=list(c("Estimated","Provided"),
                                             c("nParamInternal","nParamXreg",
                                               "nParamIntermittent","nParamAll")));

    ##### model for VES #####
    if(any(smoothType==c("ves","vets"))){
        if(!is.character(model)){
            stop(paste0("Something strange is provided instead of character object in model: ",
                        paste0(model,collapse=",")),
                 call.=FALSE);
        }

        if(length(model)==1){
            # If chosen model is "AAdN" or anything like that, we are taking the appropriate values
            if(nchar(model)==4){
                Etype <- substring(model,1,1);
                Ttype <- substring(model,2,2);
                Stype <- substring(model,4,4);
                damped <- TRUE;
                if(substring(model,3,3)!="d"){
                    warning(paste0("You have defined a strange model: ",model));
                    sowhat(model);
                    model <- paste0(Etype,Ttype,"d",Stype);
                }
                modelDo <- "estimate";
            }
            else if(nchar(model)==3){
                Etype <- substring(model,1,1);
                Ttype <- substring(model,2,2);
                Stype <- substring(model,3,3);
                damped <- FALSE;
                modelDo <- "estimate";
            }
            else{
                warning(paste0("You have defined a strange model: ",model));
                sowhat(model);
                warning("Switching to 'PPP'");
                model <- "PPP";

                Etype <- "P";
                Ttype <- "P";
                Stype <- "P";
                damped <- TRUE;
                modelDo <- "select";
            }
            modelsPool <- NULL;
        }
        else{
            modelsPool <- model;
            Etype <- "P";
            Ttype <- "P";
            Stype <- "P";
            damped <- TRUE;
            modelDo <- "select";
        }

        if((any(Etype==c("P","X","Y")) || any(Ttype==c("P","X","Y")) || any(Stype==c("P","X","Y")))){
            modelDo[] <- "select";
            if(Ttype!="N"){
                damped <- TRUE
            }
        }

        #### Check error type ####
        if(all(Etype!=c("A","M","L","P","X","Y"))){
            warning(paste0("Wrong error type: ",Etype,". Should be 'A', 'M', 'P', 'X' or 'Y'.\n",
                           "Changing to 'P'"),
                    call.=FALSE);
            Etype <- "P";
        }
        else if(Etype=="X"){
            Etype <- "A";
        }
        else if(Etype=="Y"){
            Etype <- "M";
        }

        #### Check trend type ####
        if(all(Ttype!=c("N","A","M","P","X","Y"))){
            warning(paste0("Wrong trend type: ",Ttype,". Should be 'N', 'A', 'M', 'P', 'X' or 'Y'.\n",
                           "Changing to 'P'"),
                    call.=FALSE);
            Ttype <- "P";
        }

        #### Check seasonality type ####
        # Check if seasonality makes sense
        if(all(Stype!=c("N","A","M","P","X","Y"))){
            warning(paste0("Wrong seasonality type: ",Stype,". Should be 'N', 'A', 'M', 'P', 'X' or 'Y'.",
                           "Setting to 'P'."),
                    call.=FALSE);
            if(all(lagsModel==1)){
                Stype <- "N";
            }
            else{
                Stype <- "P";
            }
        }
        if(Stype!="N" & all(lagsModel==1)){
            warning(paste0("Cannot build the seasonal model on data with frequency 1. ",
                           "Switching to non-seasonal model: ETS(",substring(model,1,nchar(model)-1),"N)"),
                    call.=FALSE);
            Stype <- "N";
        }

        if(Stype=="N"){
            initialSeason <- NULL;
            modelIsSeasonal <- FALSE;
        }
        else{
            modelIsSeasonal <- TRUE;
        }

        lagsModelMax <- max(lagsModel) * modelIsSeasonal + 1 * (!modelIsSeasonal);

        # Define the number of rows that should be in the matVt
        obsStates <- max(obsAll + lagsModelMax, obsInSample + 2*lagsModelMax);

        nComponentsNonSeasonal <- 1 + (Ttype!="N")*1;
        nComponentsAll <- nComponentsNonSeasonal + modelIsSeasonal*1;
    }

    ##### ovesModel #####
    if(is.oves(occurrence)){
        ovesModel <- occurrence$model;
        occurrence <- occurrence$occurrence;
        occurrenceModelProvided <- TRUE;
    }
    else{
        ovesModel <- model;
        occurrenceModelProvided <- FALSE;
    }

    ##### occurrence #####
    occurrence <- substring(occurrence[1],1,1);
    if(occurrence!="n"){
        ot <- (yInSample!=0)*1;
        # Matrix of non-zero observations for the loss function
        otObs <- diag(rowSums(ot));
        for(i in 1:nSeries){
            for(j in 1:nSeries){
                if(i==j){
                    next;
                }
                otObs[i,j] <- min(otObs[i,i],otObs[j,j]);
            }
        }
    }
    else{
        ot <- matrix(1,nrow=nrow(yInSample),ncol=ncol(yInSample));
        otObs <- matrix(obsInSample,nSeries,nSeries);
        ovesModel <- NULL;
    }

    # If the data is not intermittent, let's assume that the parameter was switched unintentionally.
    if(all(ot==1) & occurrence!="n"){
        occurrence <- "n";
        occurrenceModelProvided <- FALSE;
    }
    # Boolean for the occurrence
    if(occurrence=="n"){
        occurrenceModel <- FALSE;
    }
    else{
        occurrenceModel <- TRUE;
    }

    # Check if multiplicative models can be fitted
    allowMultiplicative <- !((any(yInSample<=0) && !occurrenceModel) || (occurrenceModel && any(yInSample<0)));
    modelIsMultiplicative <- FALSE;

    # Check if multiplicative model can be applied
    if(any(c(Etype,Ttype,Stype) %in% c("M","P","Y"))){
        if(allowMultiplicative){
            if(any(c(Etype,Ttype,Stype) %in% c("A","X"))){
                warning("Mixed models are not available. Switching to pure multiplicative.",call.=FALSE);
                Etype <- switch(Etype,
                                "A"="M",
                                "X"="Y",
                                Etype);
                Ttype <- switch(Ttype,
                                "A"="M",
                                "X"="Y",
                                Ttype);
                Stype <- switch(Stype,
                                "A"="M",
                                "X"="Y",
                                Stype);
            }
            # If the components are non-additive, then the model is multiplicative
            # if(all(c(Etype,Ttype,Stype)!="A") || all(c(Etype,Ttype,Stype)!="X")){
            #     yInSample <- log(yInSample);
            #     modelIsMultiplicative[] <- TRUE;
            # }
        }
        else{
            if(!occurrenceModel){
                warning("Sorry, but we cannot construct multiplicative model on non-positive data. ",
                        "Changing to additive.",
                        call.=FALSE);
                Etype <- "A";
                Ttype <- switch(Ttype,"M"="A","Y"=,"P"="X",Ttype);
                Stype <- switch(Stype,"M"="A","Y"=,"P"="X",Stype);
                modelIsMultiplicative[] <- FALSE;
            }
            else{
                # yInSample[ot==1] <- log(yInSample[ot==1]);
                Etype <- "M";
                Ttype <- ifelse(Ttype=="A","M",Ttype);
                Stype <- ifelse(Stype=="A","M",Stype);
                modelIsMultiplicative[] <- TRUE;
            }
        }
    }

    # Binaries for trend and seasonal
    modelIsTrendy <- Ttype!="N";
    modelIsSeasonal <- Stype!="N";

    # This is the number of parameters to estimate per series
    nParamMax <- 0;

    #### VES parameters ####
    if(smoothType=="ves"){
        ##### * Persistence matrix ####
        # persistence type can be: "i" - individual, "d" - dependent, "c" - common (all),
        # "s" - seasonal smoothing parameter is the same
        persistenceValue <- persistence;
        if(is.null(persistenceValue)){
            warning("persistence value is not selected. Switching to group.");
            persistenceType <- "c";
            persistenceEstimate <- TRUE;
        }
        else{
            if(is.character(persistenceValue)){
                persistenceValue <- substring(persistenceValue[1],1,1);
                if(all(persistenceValue!=c("c","i","d","s"))){
                    warning("You asked for a strange persistence value. We don't do that here. ",
                            "Switching to group",
                            call.=FALSE);
                    persistenceType <- "c";
                }
                else{
                    if(persistenceValue=="s" & Stype=="N"){
                        warning(paste0("Non-seasonal model is selected, but you've asked for common ",
                                       "seasonal smoothing parameter. ",
                                       "Switching to persistence='individual'."),
                                call.=FALSE);
                        persistenceValue <- "i";
                    }
                    persistenceType <- persistenceValue;
                }
                persistenceValue <- NULL;
                persistenceEstimate <- TRUE;
            }
            else if(is.numeric(persistenceValue)){
                if(all(length(persistenceValue) != c(nComponentsAll*nSeries^2,nComponentsAll))){
                    warning(paste0("Length of persistence matrix is wrong! It should be either ",
                                   nComponentsAll*nSeries^2, " or ", nComponentsAll,
                                   " instead of ",length(persistenceValue),".\n",
                                   "Values of persistence matrix will be estimated as group."),call.=FALSE);
                    persistenceValue <- NULL;
                    persistenceType <- "c";
                    persistenceEstimate <- TRUE;
                }
                else{
                    persistenceType <- "p";
                    persistenceEstimate <- FALSE;

                    if(length(persistenceValue)==nComponentsAll){
                        persistenceBuffer <- matrix(0,nSeries*nComponentsAll,nSeries);
                        for(i in 1:nSeries){
                            persistenceBuffer[1:nComponentsAll+nComponentsAll*(i-1),i] <- persistenceValue;
                        }
                        persistenceValue <- persistenceBuffer;
                        parametersNumber[2,1] <- parametersNumber[2,1] + length(as.vector(persistenceValue));
                    }
                    else{
                        ### Check the persistence matrix in order to decide number of parameters
                        persistencePartial <- matrix(persistenceValue[1:nComponentsAll,1:nSeries],
                                                     nComponentsAll,nSeries);
                        persistenceValue <- matrix(persistenceValue,nSeries*nComponentsAll,nSeries);

                        # Check if persistence is dependent
                        if(all(persistencePartial[,nSeries]==0)){
                            # Check if persistence is grouped
                            if(persistenceValue[1,1]==persistenceValue[1+nComponentsAll,nSeries]){
                                parametersNumber[2,1] <- parametersNumber[2,1] + nSeries;
                            }
                            else{
                                parametersNumber[2,1] <- parametersNumber[2,1] + nSeries*nComponentsAll;
                            }
                        }
                        else{
                            parametersNumber[2,1] <- parametersNumber[2,1] + length(as.vector(persistenceValue));
                        }
                    }
                }
            }
            else if(!is.numeric(persistenceValue)){
                warning(paste0("persistence matrix is not numeric!\n",
                               "Values of persistence matrix will be estimated as group."),call.=FALSE);
                persistenceValue <- NULL;
                persistenceType <- "c";
                persistenceEstimate <- TRUE;
            }
        }

        # If it is individual, then it increases by nComponentsAll
        if(persistenceType=="i"){
            # if(any(persistenceType==c("c","i","s"))){
            nParamMax <- nParamMax + nComponentsAll;
        }
        # The seasonal is shared across series, the other parameters are individual
        else if(persistenceType=="s"){
            nParamMax <- nParamMax + nComponentsNonSeasonal + 1/nSeries;
        }
        # All parameters are shared
        else if(persistenceType=="c"){
            nParamMax <- nParamMax + nComponentsAll/nSeries;
        }
        else if(persistenceType=="d"){
            # In case with "dependent" the whole matrix needs to be estimated
            nParamMax <- nParamMax + nComponentsAll*nSeries;
        }

        ##### * Transition matrix ####
        # transition type can be: "i" - individual, "d" - dependent, "c" - common
        transitionValue <- transition;
        if(is.null(transitionValue)){
            warning("transition value is not selected. Switching to common");
            transitionType <- "c";
            transitionEstimate <- FALSE;
        }
        else{
            if(is.character(transitionValue)){
                transitionValue <- substring(transitionValue[1],1,1);
                if(all(transitionValue!=c("c","i","d"))){
                    warning("You asked for a strange transition value. We don't do that here. ",
                            "Switching to common",
                            call.=FALSE);
                    transitionType <- "c";
                }
                else{
                    transitionType <- transitionValue;
                }
                transitionValue <- NULL;
                transitionEstimate <- FALSE;
            }
            else if(is.numeric(transitionValue)){
                if(all(length(transitionValue) != c((nSeries*nComponentsAll)^2,nComponentsAll^2))){
                    warning(paste0("Length of transition matrix is wrong! It should be either ",
                                   (nSeries*nComponentsAll)^2, " or ", nComponentsAll^2,
                                   " instead of ",length(transitionValue),".\n",
                                   "Values of transition matrix will be estimated as a common one."),
                            call.=FALSE);
                    transitionValue <- NULL;
                    transitionType <- "c";
                    transitionEstimate <- FALSE;
                }
                else{
                    transitionType <- "p";
                    transitionEstimate <- FALSE;
                    ### Check the transition matrix in order to decide number of parameters
                    transitionPartial <- matrix(transitionValue[1:nComponentsAll,1:nComponentsAll],
                                                nComponentsAll,nComponentsAll);
                    transitionIsStandard <- FALSE;
                    transitionContainsPhi <- FALSE;
                    if(ncol(transitionPartial)==3){
                        if(all(transitionPartial[,c(1,3)]==matrix(c(1,0,0,1,1,0,0,0,1),3,3)[,c(1,3)])){
                            transitionIsStandard <- TRUE;
                            if(transitionPartial[2,2]!=1){
                                # If there is phi in the matrix, add it
                                transitionContainsPhi <- TRUE;
                            }
                        }
                    }
                    else if(ncol(transitionPartial)==2){
                        if(all(transitionPartial[,1]==c(1,0))){
                            transitionIsStandard <- TRUE;
                            if(transitionPartial[2,2]!=1){
                                # If there is phi in the matrix, add it
                                transitionContainsPhi <- TRUE;
                            }
                        }
                    }
                    else{
                        if(transitionPartial[1,1]==1){
                            transitionIsStandard <- TRUE;
                        }
                    }
                    # If transition is not standard, take unique values from it.
                    if(!transitionIsStandard){
                        parametersNumber[2,1] <- parametersNumber[2,1] + length(as.vector(transitionValue));
                    }
                    else{
                        # If there is phi, check if it is grouped
                        if(transitionContainsPhi){
                            # If phi is grouped, add one parameter
                            if(transitionValue[2,2]==transitionValue[2+nComponentsAll,2+nComponentsAll]){
                                parametersNumber[2,1] <- parametersNumber[2,1] + 1;
                                phi <- transitionValue[2,2];
                            }
                            # Else phi is individual
                            else{
                                parametersNumber[2,1] <- parametersNumber[2,1] + nSeries;
                                phi <- rep(NA,nSeries);
                                for(i in 1:nSeries){
                                    phi[i] <- transitionValue[2+(i-1)*nComponentsAll,2+(i-1)*nComponentsAll];
                                }
                            }
                        }
                    }

                    if(length(transitionValue) == nComponentsAll^2){
                        transitionValue <- matrix(transitionValue,nComponentsAll,nComponentsAll);
                        transitionBuffer <- diag(nSeries*nComponentsAll);
                        for(i in 1:nSeries){
                            transitionBuffer[c(1:nComponentsAll)+nComponentsAll*(i-1),
                                             c(1:nComponentsAll)+nComponentsAll*(i-1)] <- transitionValue;
                        }
                        transitionValue <- transitionBuffer;
                    }
                    else{
                        transitionValue <- matrix(transitionValue,nSeries*nComponentsAll,nSeries*nComponentsAll);
                    }
                }
            }
            else if(!is.numeric(transitionValue)){
                warning(paste0("transition matrix is not numeric!\n",
                               "Values of transition vector will be estimated as a common group."),
                        call.=FALSE);
                transitionValue <- NULL;
                transitionType <- "c";
                transitionEstimate <- FALSE;
            }
        }

        if(transitionType=="d"){
            ## !!! Each separate transition matrix is not evaluated, but the off-diagonals are
            transitionEstimate <- TRUE;
            nParamMax <- nParamMax + (nSeries-1)*nSeries*nComponentsAll^2;
        }

        ##### * Damping parameter ####
        # phi type can be: "i" - individual, "c" - common
        dampedValue <- phi;
        if(transitionType!="p"){
            if(damped){
                if(is.null(dampedValue)){
                    warning("phi value is not selected. Switching to common.");
                    dampedType <- "c";
                    dampedEstimate <- TRUE;
                }
                else{
                    if(is.character(dampedValue)){
                        dampedValue <- substring(dampedValue[1],1,1);
                        if(all(dampedValue!=c("i","c"))){
                            warning("You asked for a strange phi value. We don't do that here. ",
                                    "Switching to common.",
                                    call.=FALSE);
                            dampedType <- "c";
                        }
                        else{
                            dampedType <- dampedValue;
                        }
                        dampedValue <- matrix(1,nSeries,1);
                        dampedEstimate <- TRUE;
                    }
                    else if(is.numeric(dampedValue)){
                        if((length(dampedValue) != nSeries) & (length(dampedValue)!= 1)){
                            warning(paste0("Length of phi vector is wrong! It should be ",
                                           nSeries,
                                           " instead of ",length(dampedValue),".\n",
                                           "Values of phi vector will be estimated as a common one."),
                                    call.=FALSE);
                            dampedValue <- matrix(1,nSeries,1);
                            dampedType <- "c";
                            dampedEstimate <- TRUE;
                        }
                        else{
                            dampedType <- "p";
                            dampedValue <- matrix(dampedValue,nSeries,1);
                            dampedEstimate <- FALSE;
                            parametersNumber[2,1] <- parametersNumber[2,1] + length(as.vector(dampedValue));
                        }
                    }
                    else if(!is.numeric(dampedValue)){
                        warning(paste0("phi vector is not numeric!\n",
                                       "Values of phi vector will be estimated as a common one."),
                                call.=FALSE);
                        dampedValue <- matrix(1,nSeries,1);
                        dampedType <- "c";
                        dampedEstimate <- TRUE;
                    }
                }

                if(any(dampedType==c("c","i"))){
                    dampedValue <- matrix(1,nSeries,1);
                    # In case of common, the parameter is shared.
                    if(dampedType=="c"){
                        nParamMax <- nParamMax + 1/nSeries;
                    }
                    else{
                        nParamMax <- nParamMax + 1;
                    }
                }
            }
            else{
                dampedValue <- matrix(1,nSeries,1);
                dampedType <- "c";
                dampedEstimate <- FALSE;
            }
        }
        else{
            dampedType <- "c";
            dampedEstimate <- FALSE;
        }

        ##### * initials ####
        # initial type can be: "i" - individual, "c" - common
        initialValue <- initial;
        if(is.null(initialValue)){
            warning("Initial value is not selected. Switching to individual.");
            initialType <- "i";
            initialEstimate <- TRUE;
        }
        else{
            if(is.character(initialValue)){
                initialValue <- substring(initialValue[1],1,1);
                if(all(initialValue!=c("i","c"))){
                    warning("You asked for a strange initial value. We don't do that here. ",
                            "Switching to individual.",
                            call.=FALSE);
                    initialType <- "i";
                }
                else{
                    initialType <- initialValue;
                }
                initialValue <- NULL;
                initialEstimate <- TRUE;
            }
            else if(is.numeric(initialValue)){
                if(length(initialValue)>2*nSeries){
                    warning(paste0("Length of initial vector is wrong! It should not be greater than",
                                   2*nSeries,"\n",
                                   "Values of initial vector will be estimated."),call.=FALSE);
                    initialValue <- NULL;
                    initialType <- "i";
                    initialEstimate <- TRUE;
                }
                else{
                    if(all(length(initialValue) != c(nComponentsNonSeasonal,nComponentsNonSeasonal * nSeries))){
                        warning(paste0("Length of initial vector is wrong! It should be either ",
                                       nComponentsNonSeasonal*nSeries, " or ", nComponentsNonSeasonal,
                                       " instead of ",length(initialValue),".\n",
                                       "Values of initial vector will be estimated."),call.=FALSE);
                        initialValue <- NULL;
                        initialType <- "i";
                        initialEstimate <- TRUE;
                    }
                    else{
                        initialType <- "p";
                        initialValue <- matrix(initialValue,nComponentsNonSeasonal * nSeries,1);
                        initialEstimate <- FALSE;
                        parametersNumber[2,1] <- parametersNumber[2,1] + length(as.vector(initialValue));
                    }
                }
            }
            else if(!is.numeric(initialValue)){
                warning(paste0("Initial vector is not numeric!\n",
                               "Values of initial vector will be estimated."),call.=FALSE);
                initialValue <- NULL;
                initialType <- "i";
                initialEstimate <- TRUE;
            }
        }

        # Individual initials
        if(initialType=="i"){
            nParamMax <- nParamMax + nComponentsNonSeasonal;
        }
        # Common initials are shared across series
        else{
            nParamMax <- nParamMax + nComponentsNonSeasonal / nSeries;
        }

        ##### * initialSeason and seasonal for VES #####
        # Here we should check if initialSeason is character or not...
        # if length(initialSeason) == yFrequency*nSeries, then ok
        # if length(initialSeason) == yFrequency, then use it for all nSeries
        if(Stype!="N"){
            #### Seasonal component
            seasonalType <- seasonal;
            if(is.null(seasonalType)){
                warning("The type of the seasonal component is not selected. Switching to individual.");
                seasonalType <- "i";
            }
            else{
                if(is.character(seasonalType)){
                    seasonalType <- substring(seasonalType[1],1,1);
                    if(all(seasonalType!=c("i","c"))){
                        warning("You asked for a strange seasonal value. We don't do that here. ",
                                "Switching to individual.",
                                call.=FALSE);
                        seasonalType <- "i";
                    }
                    else if(seasonalType=="c"){
                        if(Stype=="N"){
                            warning("Common seasonal model does not make sense with the non-seasonal ",
                                    "ETS. Changing to individual",
                                    call.=FALSE);
                            seasonalType <- "i";
                        }
                        else{
                            # If the transition matrix provide is full, cut off all the seasonals
                            # except for the first one.
                            if(transitionType=="p" && nrow(transitionValue)==nSeries*nComponentsAll){
                                warning(paste0("The transition matrix you provided contains too many rows ",
                                               "for the common seasonal model.",
                                               "Using only the first seasonal one."), call.=FALSE);
                                transitionValue <- rbind(cbind(transitionValue[-(c(1:nSeries)*nComponentsAll),
                                                                               -(c(1:nSeries)*nComponentsAll)],
                                                               0),
                                                         c(transitionValue[nComponentsAll*nSeries,
                                                                           -(c(1:nSeries)*nComponentsAll)],1));
                            }
                            # Do similar stuff for the persistence
                            if(persistenceType=="p" && nrow(persistenceValue)==nSeries*nComponentsAll){
                                warning(paste0("The persistence matrix you provided contains too many rows ",
                                               "for the common seasonal model.",
                                               "Using only the first seasonal one."), call.=FALSE);
                                persistenceValue <- rbind(persistenceValue[-(c(1:nSeries)*nComponentsAll),],
                                                          persistenceValue[nComponentsAll,]);
                            }
                        }
                    }
                }
                else{
                    warning("A weird stuff is provided for the seasonal component. Switching to individual.",
                            call.=FALSE);
                    seasonalType <- "i";
                }
            }

            #### Initials
            initialSeasonValue <- initialSeason;
            if(is.null(initialSeasonValue)){
                warning("Initial value is not selected. Switching to common");
                initialSeasonType <- "c";
                initialSeasonEstimate <- TRUE;
            }
            else{
                if(is.character(initialSeasonValue)){
                    initialSeasonValue <- substring(initialSeasonValue[1],1,1);
                    if(all(initialSeasonValue!=c("i","c"))){
                        warning("You asked for a strange initialSeason value. We don't do that here. ",
                                "Switching to common",
                                call.=FALSE);
                        initialSeasonType <- "c";
                    }
                    else{
                        if(seasonalType=="c" && initialSeasonValue=="i"){
                            warning("initialSeason='i' does not work with seasonalType='c'. ",
                                    "Switching to common.",
                                    call.=FALSE);
                            initialSeasonValue <- "c";
                        }
                        initialSeasonType <- initialSeasonValue;
                    }
                    initialSeasonValue <- NULL;
                    initialSeasonEstimate <- TRUE;
                }
                else if(is.numeric(initialSeasonValue)){
                    if(all(length(initialSeasonValue)!=c(lagsModelMax,lagsModelMax*nSeries))){
                        warning(paste0("The length of initialSeason is wrong! ",
                                       "It should correspond to the frequency of the data.",
                                       "Values of initialSeason will be estimated as a common one."),
                                call.=FALSE);
                        initialSeasonValue <- NULL;
                        initialSeasonType <- "c";
                        initialSeasonEstimate <- TRUE;
                    }
                    else{
                        if(seasonalType=="i"){
                            initialSeasonValue <- matrix(initialSeasonValue,nSeries,lagsModelMax);
                        }
                        else{
                            if(length(initialSeasonValue)!=lagsModelMax){
                                warning(paste0("The initialSeason you provided contains too many elements ",
                                               "for the common seasonal model.",
                                               "Using only the first ",lagsModelMax," values."), call.=FALSE);
                            }
                            initialSeasonValue <- matrix(initialSeasonValue[1:lagsModelMax],1,lagsModelMax);
                        }
                        initialSeasonType <- "p";
                        initialSeasonEstimate <- FALSE;
                        parametersNumber[2,1] <- parametersNumber[2,1] + length(as.vector(initialSeasonValue));
                    }
                }
                else if(!is.numeric(initialSeasonValue)){
                    warning(paste0("Initial vector is not numeric!\n",
                                   "Values of initialSeason vector will be estimated as a common one."),
                            call.=FALSE);
                    initialSeasonValue <- NULL;
                    initialSeasonType <- "c";
                    initialSeasonEstimate <- TRUE;
                }
            }

            if(initialSeasonType=="i"){
                nParamMax <- nParamMax + lagsModelMax;
            }
            else if(initialSeasonType=="c"){
                nParamMax <- nParamMax + lagsModelMax / nSeries;
            }
        }
        else{
            initialSeasonValue <- NULL;
            initialSeasonType <- "c";
            initialSeasonEstimate <- FALSE;
            seasonalType <- "i";
        }
    }
    #### VETS parameters ####
    # vets doesn't have provided parameters.
    else if(smoothType=="vets"){
        # Define common parameters
        parameters <- substr(parameters,1,1);
        # Define common initials
        initials <- substr(initials,1,1);
        # Defin common components
        components <- substr(components,1,1);
    }

    ##### Loss function type #####
    loss <- match.arg(loss, c("likelihood","diagonal","trace"));

    # Modify loss for the oves model
    if(Etype=="L"){
        loss[] <- "occurrence";
    }

    # If it is likelihood, we also need to estimate the full covariance matrix
    if(loss=="likelihood"){
        nParamMax <- nParamMax + nSeries;
    }
    # Otherwise, these are just variances of the data
    else{
        nParamMax <- nParamMax + 1;
    }

    normalizer <- sum(colMeans(abs(diff(t(yInSample))),na.rm=TRUE));

    ##### Information Criteria #####
    ic <- ic[1];
    if(all(ic!=c("AICc","AIC","BIC","BICc"))){
        warning(paste0("Strange type of information criteria defined: ",ic,". Switching to 'AICc'."),
                call.=FALSE);
        ic <- "AICc";
    }

    ##### interval, intervalType, level #####
    # intervalType <- interval[1];
    # Check the provided type of interval

    # if(is.logical(intervalType)){
    #     if(intervalType){
    #         intervalType <- "i";
    #     }
    #     else{
    #         intervalType <- "none";
    #     }
    # }

    # if(all(intervalType!=c("c","u","i","l","n","none","conditional","unconditional",
    #                        "individual","likelihood"))){
    #     warning(paste0("Wrong type of interval: '",intervalType, "'. Switching to 'conditional'."),
    #             call.=FALSE);
    #     intervalType <- "i";
    # }
    #
    # if(intervalType=="none"){
    #     intervalType <- "n";
    #     interval <- FALSE;
    # }
    # else if(intervalType=="conditional"){
    #     intervalType <- "c";
    #     interval <- TRUE;
    # }
    # else if(intervalType=="unconditional"){
    #     intervalType <- "u";
    #     interval <- TRUE;
    # }
    # else if(intervalType=="individual"){
    #     intervalType <- "i";
    #     interval <- TRUE;
    # }
    # else if(intervalType=="likelihood"){
    #     intervalType <- "l";
    #     interval <- TRUE;
    # }
    # else{
    #     interval <- TRUE;
    # }
    #
    # if(level>1){
    #     level <- level / 100;
    # }

    ##### bounds #####
    bounds <- substring(bounds[1],1,1);
    if(all(bounds!=c("u","a","n"))){
        warning("Strange bounds are defined. Switching to 'admissible'.",call.=FALSE);
        bounds <- "a";
    }

    ##### Check number of observations vs number of max parameters #####
    if(obsInSample <= nParamMax){
        stop(paste0("Not enough observations for the reasonable fit. Number of parameters per series is ",
                    nParamMax," while the number of observations is ",obsInSample,"."),call.=FALSE);
    }

    ##### Fisher Information #####
    if(!exists("FI")){
        FI <- FALSE;
    }
    else{
        if(!is.logical(FI)){
            FI <- FALSE;
        }
        if(!requireNamespace("numDeriv",quietly=TRUE) & FI){
            warning(paste0("Sorry, but you don't have 'numDeriv' package, ",
                           "which is required in order to produce Fisher Information.",call.=FALSE));
            FI <- FALSE;
        }
    }

    ##### Ellipsis thingies #####
    if(!is.null(ellipsis$B)){
        B <- ellipsis$B;
    }
    else{
        B <- NULL;
    }
    if(!is.null(ellipsis$ub)){
        ub <- ellipsis$ub;
    }
    else{
        ub <- NULL;
    }
    if(!is.null(ellipsis$lb)){
        lb <- ellipsis$lb;
    }
    else{
        lb <- NULL;
    }
    if(!is.null(ellipsis$maxeval)){
        maxeval <- ellipsis$maxeval;
    }
    else{
        maxeval <- NULL;
    }
    if(!is.null(ellipsis$algorithm1)){
        algorithm1 <- ellipsis$algorithm1;
    }
    else{
        algorithm1 <- "NLOPT_LN_BOBYQA";
    }
    if(!is.null(ellipsis$algorithm2)){
        algorithm2 <- ellipsis$algorithm2;
    }
    else{
        algorithm2 <- "NLOPT_LN_NELDERMEAD";
    }
    if(!is.null(ellipsis$xtol_rel1)){
        xtol_rel1 <- ellipsis$xtol_rel1;
    }
    else{
        xtol_rel1 <- 1e-8;
    }
    if(!is.null(ellipsis$xtol_rel2)){
        xtol_rel2 <- ellipsis$xtol_rel2;
    }
    else{
        xtol_rel2 <- 1e-6;
    }
    if(!is.null(ellipsis$print_level)){
        print_level <- ellipsis$print_level;
    }
    else{
        print_level <- 0;
    }

    ##### Return values to previous environment #####
    assign("h",h,ParentEnvironment);
    assign("silent",silent,ParentEnvironment);
    assign("obsInSample",obsInSample,ParentEnvironment);
    assign("obsAll",obsAll,ParentEnvironment);
    assign("obsStates",obsStates,ParentEnvironment);
    assign("nSeries",nSeries,ParentEnvironment);
    assign("nParamMax",nParamMax,ParentEnvironment);
    assign("data",data,ParentEnvironment);
    assign("yInSample",yInSample,ParentEnvironment);

    # ts / zoo elements
    assign("yClasses",yClasses,ParentEnvironment);
    assign("yIndex",yIndex,ParentEnvironment);
    assign("yInSampleIndex",yInSampleIndex,ParentEnvironment);
    assign("yForecastIndex",yForecastIndex,ParentEnvironment);
    assign("yIndexAll",yIndexAll,ParentEnvironment);
    assign("yFrequency",yFrequency,ParentEnvironment);
    assign("yStart",yStart,ParentEnvironment);
    assign("yForecastStart",yForecastStart,ParentEnvironment);
    assign("yDeltat",yDeltat,ParentEnvironment);

    assign("dataNames",dataNames,ParentEnvironment);
    assign("parametersNumber",parametersNumber,ParentEnvironment);

    assign("model",model,ParentEnvironment);
    assign("modelsPool",modelsPool,ParentEnvironment);
    assign("Etype",Etype,ParentEnvironment);
    assign("Ttype",Ttype,ParentEnvironment);
    assign("Stype",Stype,ParentEnvironment);
    assign("lagsModelMax",lagsModelMax,ParentEnvironment);
    assign("modelIsTrendy",modelIsTrendy,ParentEnvironment);
    assign("modelIsSeasonal",modelIsSeasonal,ParentEnvironment);
    assign("allowMultiplicative",allowMultiplicative,ParentEnvironment);
    assign("modelIsMultiplicative",modelIsMultiplicative,ParentEnvironment);
    assign("nComponentsAll",nComponentsAll,ParentEnvironment);
    assign("nComponentsNonSeasonal",nComponentsNonSeasonal,ParentEnvironment);
    assign("damped",damped,ParentEnvironment);
    assign("modelDo",modelDo,ParentEnvironment);

    if(smoothType=="ves"){
        assign("persistenceValue",persistenceValue,ParentEnvironment);
        assign("persistenceType",persistenceType,ParentEnvironment);
        assign("persistenceEstimate",persistenceEstimate,ParentEnvironment);

        assign("transitionValue",transitionValue,ParentEnvironment);
        assign("transitionType",transitionType,ParentEnvironment);
        assign("transitionEstimate",transitionEstimate,ParentEnvironment);

        assign("dampedValue",dampedValue,ParentEnvironment);
        assign("dampedType",dampedType,ParentEnvironment);
        assign("dampedEstimate",dampedEstimate,ParentEnvironment);

        assign("initialValue",initialValue,ParentEnvironment);
        assign("initialType",initialType,ParentEnvironment);
        assign("initialEstimate",initialEstimate,ParentEnvironment);

        assign("initialSeasonValue",initialSeasonValue,ParentEnvironment);
        assign("initialSeasonType",initialSeasonType,ParentEnvironment);
        assign("initialSeasonEstimate",initialSeasonEstimate,ParentEnvironment);

        assign("seasonalType",seasonalType,ParentEnvironment);
    }
    else{
        assign("parameters",parameters,ParentEnvironment);
        assign("initials",initials,ParentEnvironment);
        assign("components",components,ParentEnvironment);
    }

    assign("loss",loss,ParentEnvironment);
    assign("normalizer",normalizer,ParentEnvironment);

    assign("ic",ic,ParentEnvironment);

    assign("occurrence",occurrence,ParentEnvironment);
    assign("ot",ot,ParentEnvironment);
    assign("otObs",otObs,ParentEnvironment);
    assign("ovesModel",ovesModel,ParentEnvironment);
    assign("occurrenceModelProvided",occurrenceModelProvided,ParentEnvironment);

    # assign("yot",yot,ParentEnvironment);
    # assign("pt",pt,ParentEnvironment);
    # assign("pt.for",pt.for,ParentEnvironment);
    # assign("nParamIntermittent",nParamIntermittent,ParentEnvironment);
    # assign("iprob",iprob,ParentEnvironment);

    assign("bounds",bounds,ParentEnvironment);

    # Stuff in ellipsis
    assign("FI",FI,ParentEnvironment);
    assign("B",B,ParentEnvironment);
    assign("ub",ub,ParentEnvironment);
    assign("lb",lb,ParentEnvironment);
    assign("maxeval",maxeval,ParentEnvironment);
    assign("algorithm1",algorithm1,ParentEnvironment);
    assign("algorithm2",algorithm2,ParentEnvironment);
    assign("xtol_rel1",xtol_rel1,ParentEnvironment);
    assign("xtol_rel2",xtol_rel2,ParentEnvironment);
    assign("print_level",print_level,ParentEnvironment);
}

##### CF for scale calculation #####
# This is needed for the models with multiplicative error term
scalerCF <- function(A, errors, scaleValue, yInSampleSum, loss){
    # Fill in the matrix
    if(loss=="likelihood"){
        scaleValue[upper.tri(scaleValue,diag=TRUE)] <- A;
        scaleValue[lower.tri(scaleValue)] <- t(scaleValue)[lower.tri(scaleValue)];
    }
    else{
        diag(scaleValue) <- A;
    }
    # If detereminant is positive, return logLik
    if(det(scaleValue)<=0){
        return(1e+100);
    }
    else{
        return(-sum(dmvnormInternal(errors, -0.5*diag(scaleValue), scaleValue, log=TRUE))+yInSampleSum);
    }
}

##### *vssFitter function* #####
vssFitter <- function(...){
    ellipsis <- list(...);
    ParentEnvironment <- ellipsis[['ParentEnvironment']];

    fitting <- vFitterWrap(switch(Etype, "M"=log(yInSample), yInSample),
                           matVt, matF, matW, matG,
                           lagsModel, Etype, Ttype, Stype, ot);
    matVt[] <- fitting$matVt;
    yFitted[] <- fitting$yfit;
    errors[] <- fitting$errors;

    if(Etype=="M"){
        yFitted <- exp(yFitted);
        if(occurrence!="n"){
            yFitted[ot==0] <- 0;
        }
    }

    # Division by nSeries gives the df per series, which agrees with Lutkepohl (2005), p.75
    nParamPerSeries <- nParam / nSeries;
    # df <- (otObs - nParamPerSeries);
    # if(intervalType!="l" && any(otObs >= nParamPerSeries)){
        df <- otObs - nParamPerSeries;
    # }
    # else{
    #     df <- otObs;
    # }

    # Divide each element by each degree of freedom
    if(loss=="diagonal"){
        Sigma <- diag(rowSums(errors^2)) / diag(df);
    }
    else{
        Sigma <- (errors %*% t(errors)) / df;
    }
    rownames(Sigma) <- colnames(Sigma) <- dataNames;

    assign("matVt",matVt,ParentEnvironment);
    assign("yFitted",yFitted,ParentEnvironment);
    assign("errors",errors,ParentEnvironment);
    assign("Sigma",Sigma,ParentEnvironment);
    assign("df",df,ParentEnvironment);
}

##### *Forecaster of state space functions* #####
vssForecaster <- function(...){
    ellipsis <- list(...);
    ParentEnvironment <- ellipsis[['ParentEnvironment']];

    # if(any((otObs - nParamPerSeries)<=0)){
    #     df <- 0;
    # }
    # else{
    # Take the minimum df for the purposes of interval construction
        df <- min(df);
    # }

    PI <- NA;

    if(h>0){
        yForecast[] <- vForecasterWrap(matrix(matVt[,(obsInSample+1):(obsInSample+lagsModelMax)],
                                              ncol=lagsModelMax),
                                     matF, matW, nSeries, h, Etype, Ttype, Stype, lagsModel);
    }
    else{
        yForecast[] <- NA;
    }

    if(any(is.na(yFitted),all(is.na(yForecast),h>0))){
        warning("Something went wrong during the optimisation and NAs were produced!",
                call.=FALSE);
        warning("Please check the input and report this error to the maintainer if it persists.",
                call.=FALSE);
    }

    if(Etype=="M"){
        yForecast[] <- exp(yForecast);
        PI[] <- exp(PI);
    }

    if(occurrence!="n"){
        if(!occurrenceModelProvided){
            ovesModel <- oves(ts(t(ot),frequency=yFrequency),
                              occurrence=occurrence, h=h, holdout=FALSE,
                              probability="dependent", model=ovesModel);
        }
        yForecast[] <- yForecast * t(ovesModel$forecast);
    }

    assign("yForecast",yForecast,ParentEnvironment);
    assign("PI",PI,ParentEnvironment);
    assign("ovesModel",ovesModel,ParentEnvironment);
}
