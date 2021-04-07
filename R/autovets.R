utils::globalVariables(c("parameters","initials","components"))

#' @importFrom utils combn
#' @rdname vets
#' @export
auto.vets <- function(y, model="PPP", lags=c(frequency(y)),
                      loss=c("likelihood","diagonal","trace"),
                      ic=c("AICc","AIC","BIC","BICc"), h=10, holdout=FALSE,
                      interval=c("none","conditional","unconditional","individual","likelihood"), level=0.95,
                      occurrence=c("none","fixed","logistic"),
                      bounds=c("admissible","usual","none"),
                      silent=TRUE, ...){
    # Copyright (C) 2021 - Inf  Ivan Svetunkov

    # Start measuring the time of calculations
    startTime <- Sys.time();

    # The function selects the restrictions on PIC elements
    ic <- match.arg(ic);

    architectorAutoVETS <- function(model){
        # Get ETS components. This is needed for modelIsTrendy and modelsIsSeasonal
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
        }
        else if(nchar(model)==3){
            Etype <- substring(model,1,1);
            Ttype <- substring(model,2,2);
            Stype <- substring(model,3,3);
            damped <- FALSE;
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
        }

        # Binaries for trend and seasonal
        modelIsTrendy <- Ttype!="N";
        modelIsSeasonal <- Stype!="N";

        parametersToCheck <- c("l","t"[modelIsTrendy],"s"[modelIsSeasonal],"d"[modelIsTrendy]);
        initialsToCheck <- c("l","t"[modelIsTrendy],"s"[modelIsSeasonal]);
        componentsToCheck <- c("l","t"[modelIsTrendy],"s"[modelIsSeasonal]);

        return(list(modelIsTrendy=modelIsTrendy,modelIsSeasonal=modelIsSeasonal,
                    parameters=parametersToCheck,initials=initialsToCheck,components=componentsToCheck))
    }

    # Prepare the call of vets()
    vetsCall <- list(...);
    vetsCall$y <- y;
    vetsCall$model <- model;
    vetsCall$lags <- lags;
    vetsCall$loss <- loss;
    vetsCall$ic <- ic;
    vetsCall$h <- h;
    vetsCall$holdout <- holdout;
    vetsCall$interval <- interval;
    vetsCall$level <- level;
    vetsCall$occurrence <- occurrence;
    vetsCall$bounds <- bounds
    vetsCall$silent <- TRUE;
    vetsCall$parameters <- "none";
    vetsCall$initials <- "none";
    vetsCall$components <- "none";

    if(!silent){
        cat("Selecting the best unrestricted model... \n");
    }

    # Select the model for the basic
    initialModel <- do.call("vets",vetsCall);
    vetsCall$model <- modelType(initialModel);

    # Get parameters, initials and components based on the selected model
    list2env(architectorAutoVETS(vetsCall$model),environment());

    # List of all the combinations of parameters restrictions
    parametersCombinations <- choose(length(parameters),c(1:length(parameters)));
    parametersCombinationNumber <- sum(parametersCombinations);
    parametersToCheck <- vector("list",parametersCombinationNumber);
    for(i in 1:length(parametersCombinations)){
        if(i==1){
            parametersToCheck[1:parametersCombinations[1]] <- as.list(as.data.frame(combn(parameters,i),stringsAsFactors=F));
        }
        else{
            parametersToCheck[sum(parametersCombinations[1:(i-1)])+1:parametersCombinations[i]] <-
                as.list(as.data.frame(combn(parameters,i),stringsAsFactors=F));
        }
    }

    # List of all the combinations of initials restrictions
    initialsCombinations <- choose(length(initials),c(1:length(initials)));
    initialsCombinationNumber <- sum(initialsCombinations);
    initialsToCheck <- vector("list",initialsCombinationNumber);
    for(i in 1:length(initialsCombinations)){
        if(i==1){
            initialsToCheck[1:initialsCombinations[1]] <- as.list(as.data.frame(combn(initials,i),stringsAsFactors=F));
        }
        else{
            initialsToCheck[sum(initialsCombinations[1:(i-1)])+1:initialsCombinations[i]] <-
                as.list(as.data.frame(combn(initials,i),stringsAsFactors=F));
        }
    }
    # List of all the combinations of components restrictions
    componentsCombinations <- choose(length(components),c(1:length(components)));
    componentsCombinationNumber <- sum(componentsCombinations);
    componentsToCheck <- vector("list",componentsCombinationNumber);
    for(i in 1:length(componentsCombinations)){
        if(i==1){
            componentsToCheck[1:componentsCombinations[1]] <- as.list(as.data.frame(combn(components,i),stringsAsFactors=F));
        }
        else{
            componentsToCheck[sum(componentsCombinations[1:(i-1)])+1:componentsCombinations[i]] <-
                as.list(as.data.frame(combn(components,i),stringsAsFactors=F));
        }
    }

    # Drop unreasonable components restrictions
    # e.g., if there is trend, then you cannot restrict level only
    if(modelIsTrendy){
        componentsToCheck[componentsToCheck %in% "l"][[1]] <- c("l","t");
        if(modelIsSeasonal){
            componentsToCheck[componentsToCheck %in% c("l","s")][[1]] <- c("l","t","s");
        }
        componentsToCheck <- unique(componentsToCheck);
    }

    #### Start the fitting of all the models ####
    # All the options from models to check
    nToCheck <- length(parametersToCheck)+length(initialsToCheck)+length(componentsToCheck);
    # +1 is for the initial model
    vetsModels <- vector("list",nToCheck+1);
    vetsModels[[1]] <- initialModel;
    rm(initialModel);
    if(!silent){
        cat(paste0("Initial model is VETS(",modelType(vetsModels[[1]]),"), IC is: ", round(vetsModels[[1]]$ICs[ic],3),"\n"));
    }

    if(!silent){
        cat("Testing parameters restrictions... ");
    }
    # Test the models with parameters restrictions
    j <- 2;
    for(i in 1:parametersCombinationNumber){
        if(i>1){
            cat(paste0(rep("\b",nchar(round((i-1)/parametersCombinationNumber,2)*100)+1),collapse=""));
        }
        cat(round(i/parametersCombinationNumber,2)*100,"\b%");

        vetsCall$parameters <- parametersToCheck[[i]];
        vetsModels[[j]] <- do.call("vets",vetsCall);
        j[] <- j+1;
    }
    # Which of the models has the lowest IC? It is the best one so far
    vetsICsParameters <- sapply(vetsModels[0:parametersCombinationNumber+1],"[[","ICs")[ic,];
    jBest <- which.min(vetsICsParameters);
    ICBest <- vetsICsParameters[jBest];
    if(jBest==1){
        vetsCall$parameters <- "none";
    }
    else{
        vetsCall$parameters <- parametersToCheck[[jBest-1]];
    }
    if(!silent){
        cat(paste0("\nParameters restrictions model is (",paste0(vetsCall$parameters,collapse=","),
                   "), IC is: ", round(ICBest,3),"\n"));
    }

    if(!silent){
        cat("Testing initials restrictions... ");
    }
    # Test the models with initials restrictions
    for(i in 1:initialsCombinationNumber){
        if(i>1){
            cat(paste0(rep("\b",nchar(round((i-1)/initialsCombinationNumber,2)*100)+1),collapse=""));
        }
        cat(round(i/initialsCombinationNumber,2)*100,"\b%");

        vetsCall$initials <- initialsToCheck[[i]];
        vetsModels[[j]] <- do.call("vets",vetsCall);
        j[] <- j+1;
    }
    # Find the model with the lowest IC from the new ones
    vetsICsInitials <- sapply(vetsModels[parametersCombinationNumber+1+1:initialsCombinationNumber],"[[","ICs")[ic,];
    jBestInitials <- which.min(vetsICsInitials);
    if(vetsICsInitials[jBestInitials]<ICBest){
        jBest <- jBestInitials+parametersCombinationNumber+1;
        ICBest <- vetsICsInitials[jBestInitials];
        vetsCall$initials <- initialsToCheck[[jBestInitials]];
    }
    else{
        vetsCall$initials <- "none";
    }
    if(!silent){
        cat(paste0("\nInitials restrictions model is (",paste0(vetsCall$initials,collapse=","),
                   "), IC is: ", round(ICBest,3),"\n"));
    }

    if(!silent){
        cat("Testing components restrictions... ");
    }
    # Test the models with initials restrictions
    for(i in 1:componentsCombinationNumber){
        if(i>1){
            cat(paste0(rep("\b",nchar(round((i-1)/componentsCombinationNumber,2)*100)+1),collapse=""));
        }
        cat(round(i/componentsCombinationNumber,2)*100,"\b%");

        vetsCall$components <- initialsToCheck[[i]];
        vetsModels[[j]] <- do.call("vets",vetsCall);
        j[] <- j+1;
    }
    # Find the model with the lowest IC from the new ones
    vetsICsComponents <- sapply(vetsModels[initialsCombinationNumber+parametersCombinationNumber+1+
                                               1:componentsCombinationNumber],"[[","ICs")[ic,];
    jBestComponents <- which.min(vetsICsComponents);
    if(vetsICsComponents[jBestComponents]<ICBest){
        jBest <- jBestComponents+initialsCombinationNumber+parametersCombinationNumber+1;
        ICBest <- vetsICsComponents[jBestComponents];
        vetsCall$components <- initialsToCheck[[jBestComponents]];
    }
    else{
        vetsCall$components <- "none";
    }
    if(!silent){
        cat(paste0("\nComponents restrictions model is (",paste0(vetsCall$components,collapse=","),
                   "), IC is: ", round(ICBest,3),"\n"));
    }

    vetsModels[[jBest]]$timeElapsed <- Sys.time()-startTime;
    return(vetsModels[[jBest]]);
}