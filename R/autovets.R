utils::globalVariables(c("parameters","initials","components"))

#' @param parallel If TRUE, the estimation of ADAM models is done in parallel (used in \code{auto.adam} only).
#' If the number is provided (e.g. \code{parallel=41}), then the specified number of cores is set up.
#' WARNING! Packages \code{foreach} and either \code{doMC} (Linux and Mac only)
#' or \code{doParallel} are needed in order to run the function in parallel.
#' @importFrom utils combn
#' @rdname vets
#' @export
auto.vets <- function(data, model="PPP", lags=c(frequency(data)),
                      loss=c("likelihood","diagonal","trace"),
                      ic=c("AICc","AIC","BIC","BICc"), h=10, holdout=FALSE,
                      occurrence=c("none","fixed","logistic"),
                      bounds=c("admissible","usual","none"),
                      silent=TRUE, parallel=FALSE, ...){
    # The function selects the restrictions on PIC elements
    # Copyright (C) 2021 - Inf  Ivan Svetunkov

    # Start measuring the time of calculations
    startTime <- Sys.time();

    ic <- match.arg(ic);
    loss <- match.arg(loss);

    #### Architector for ETS part ####
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

    #### Call the basic vets() ####
    # Prepare the call of vets()
    vetsCall <- list(...);
    vetsCall$data <- data;
    vetsCall$model <- model;
    vetsCall$lags <- lags;
    vetsCall$ic <- ic;
    vetsCall$h <- h;
    vetsCall$holdout <- holdout;
    vetsCall$occurrence <- occurrence;
    vetsCall$bounds <- bounds
    vetsCall$silent <- TRUE;
    vetsCall$parameters <- "none";
    vetsCall$initials <- "none";
    vetsCall$components <- "none";
    vetsCall$loss <- loss;
    # Use diagonal in order to be able to select something
    if(loss=="likelihood"){
        vetsCall$loss <- "diagonal";
    }

    if(!silent){
        cat("Selecting the best unrestricted model... \n");
    }

    # Select the model for the basic
    initialModel <- do.call("vets",vetsCall);
    vetsCall$model[] <- modelType(initialModel);
    # Use the specified loss
    vetsCall$loss[] <- loss;

    # Get parameters, initials and components based on the selected model
    list2env(architectorAutoVETS(vetsCall$model),environment());

    #### Prepare all pools ####
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
        componentsCombinationNumber <- length(componentsToCheck);
    }

    #### Parallel calculations ####
    # Check the parallel parameter and set the number of cores
    if(is.numeric(parallel)){
        nCores <- parallel;
        parallel <- TRUE
    }
    else{
        if(parallel){
            nCores <- min(parallel::detectCores() - 1, parametersCombinationNumber,
                          initialsCombinationNumber, componentsCombinationNumber);
        }
    }

    # If this is parallel, then load the required packages
    if(parallel){
        if(!requireNamespace("foreach", quietly = TRUE)){
            stop("In order to run the function in parallel, 'foreach' package must be installed.", call. = FALSE);
        }
        if(!requireNamespace("parallel", quietly = TRUE)){
            stop("In order to run the function in parallel, 'parallel' package must be installed.", call. = FALSE);
        }

        # Check the system and choose the package to use
        if(Sys.info()['sysname']=="Windows"){
            if(requireNamespace("doParallel", quietly = TRUE)){
                cat("Setting up", nCores, "clusters using 'doParallel'...\n");
                cluster <- parallel::makeCluster(nCores);
                doParallel::registerDoParallel(cluster);
            }
            else{
                stop("Sorry, but in order to run the function in parallel, you need 'doParallel' package.",
                     call. = FALSE);
            }
        }
        else{
            if(requireNamespace("doMC", quietly = TRUE)){
                doMC::registerDoMC(nCores);
                cluster <- NULL;
            }
            else if(requireNamespace("doParallel", quietly = TRUE)){
                cat("Setting up", nCores, "clusters using 'doParallel'...\n");
                cluster <- parallel::makeCluster(nCores);
                doParallel::registerDoParallel(cluster);
            }
            else{
                stop(paste0("Sorry, but in order to run the function in parallel, you need either ",
                            "'doMC' (prefered) or 'doParallel' package."),
                     call. = FALSE);
            }
        }
    }
    else{
        cluster <- NULL;
    }

    #### Start fitting of all the models ####
    # All the options from models to check
    nToCheck <- length(parametersToCheck)+length(initialsToCheck)+length(componentsToCheck);
    # +1 is for the initial model
    vetsModels <- vector("list",nToCheck+1);
    vetsModels[[1]] <- initialModel;
    rm(initialModel);
    if(!silent){
        cat(paste0("Initial model is VETS(",modelType(vetsModels[[1]]),"), ",ic," is: ", round(vetsModels[[1]]$ICs[ic],3),"\n"));
    }

    if(!silent){
        cat("Testing initials restrictions... ");
    }
    if(!parallel){
        j <- 2;
        # Test the models with initials restrictions
        for(i in 1:initialsCombinationNumber){
            if(!silent){
                if(i>1){
                    cat(paste0(rep("\b",nchar(round((i-1)/initialsCombinationNumber,2)*100)+1),collapse=""));
                }
                cat(round(i/initialsCombinationNumber,2)*100,"\b%");
            }

            vetsCall$initials <- initialsToCheck[[i]];
            vetsModels[[j]] <- do.call("vets",vetsCall);
            j[] <- j+1;
        }
    }
    else{
        vetsModels[1+1:initialsCombinationNumber] <-
            foreach::`%dopar%`(foreach::foreach(i=1:initialsCombinationNumber),{
            vetsCall$initials <- initialsToCheck[[i]];
            return(do.call("vets",vetsCall));
        })
    }

    # Which of the models has the lowest IC? It is the best one so far
    vetsICsInitials <- sapply(vetsModels[1+0:initialsCombinationNumber],"[[","ICs")[ic,];
    jBest <- which.min(vetsICsInitials);
    ICBest <- vetsICsInitials[jBest];
    if(jBest==1){
        vetsCall$initials <- "none";
    }
    else{
        vetsCall$initials <- initialsToCheck[[jBest-1]];
    }
    if(!silent){
        cat(paste0("\nInitials restrictions model is (",paste0(vetsCall$initials,collapse=","),
                   "), ",ic," is: ", round(ICBest,3),"\n"));
    }

    if(!silent){
        cat("Testing parameters restrictions... ");
    }
    if(!parallel){
        # Test the models with parameters restrictions
        for(i in 1:parametersCombinationNumber){
            if(!silent){
                if(i>1){
                    cat(paste0(rep("\b",nchar(round((i-1)/parametersCombinationNumber,2)*100)+1),collapse=""));
                }
                cat(round(i/parametersCombinationNumber,2)*100,"\b%");
            }

            vetsCall$parameters <- parametersToCheck[[i]];
            vetsModels[[j]] <- do.call("vets",vetsCall);
            j[] <- j+1;
        }
    }
    else{
        vetsModels[1+initialsCombinationNumber+1:parametersCombinationNumber] <-
            foreach::`%dopar%`(foreach::foreach(i=1:parametersCombinationNumber),{
            vetsCall$parameters <- parametersToCheck[[i]];
            return(do.call("vets",vetsCall));
        })
    }
    # Find the model with the lowest IC from the new ones
    vetsICsParameters <- sapply(vetsModels[1+initialsCombinationNumber+1:parametersCombinationNumber],"[[","ICs")[ic,];
    jBestParameters <- which.min(vetsICsParameters);
    if(vetsICsParameters[jBestParameters]<ICBest){
        jBest <- jBestParameters+initialsCombinationNumber+1;
        ICBest <- vetsICsParameters[jBestParameters];
        vetsCall$parameters <- parametersToCheck[[jBestParameters]];
    }
    else{
        vetsCall$parameters <- "none";
    }
    if(!silent){
        cat(paste0("\nParameters restrictions model is (",paste0(vetsCall$parameters,collapse=","),
                   "), ",ic," is: ", round(ICBest,3),"\n"));
    }

    if(!silent){
        cat("Testing components restrictions... ");
    }
    if(!parallel){
        # Test the models with initials restrictions
        for(i in 1:componentsCombinationNumber){
            if(!silent){
                if(i>1){
                    cat(paste0(rep("\b",nchar(round((i-1)/componentsCombinationNumber,2)*100)+1),collapse=""));
                }
                cat(round(i/componentsCombinationNumber,2)*100,"\b%");
            }

            vetsCall$components <- componentsToCheck[[i]];
            vetsModels[[j]] <- do.call("vets",vetsCall);
            j[] <- j+1;
        }
    }
    else{
        vetsModels[1+initialsCombinationNumber+parametersCombinationNumber+1:componentsCombinationNumber] <-
            foreach::`%dopar%`(foreach::foreach(i=1:componentsCombinationNumber),{
            vetsCall$components <- componentsToCheck[[i]];
            return(do.call("vets",vetsCall));
        })
    }
    # Find the model with the lowest IC from the new ones
    vetsICsComponents <- sapply(vetsModels[1+initialsCombinationNumber+parametersCombinationNumber+
                                               1:componentsCombinationNumber],"[[","ICs")[ic,];
    jBestComponents <- which.min(vetsICsComponents);
    if(vetsICsComponents[jBestComponents]<ICBest){
        jBest <- jBestComponents+1+initialsCombinationNumber+parametersCombinationNumber;
        ICBest <- vetsICsComponents[jBestComponents];
        vetsCall$components <- componentsToCheck[[jBestComponents]];
    }
    else{
        vetsCall$components <- "none";
    }
    if(!silent){
        cat(paste0("\nComponents restrictions model is (",paste0(vetsCall$components,collapse=","),
                   "), ",ic," is: ", round(ICBest,3),"\n"));
    }

    # Check if the clusters have been made
    if(!is.null(cluster)){
        parallel::stopCluster(cluster);
    }

    vetsModels[[jBest]]$timeElapsed <- Sys.time()-startTime;

    return(vetsModels[[jBest]]);
}