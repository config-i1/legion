#' @importFrom stats logLik
#' @export
logLik.legion <- function(object,...){
    obs <- nobs(object);
    nParamPerSeries <- nparam(object);
    structure(object$logLik,nobs=obs,df=nParamPerSeries,class="logLik");
}

#' @export
logLik.oves <- function(object,...){
    obs <- nobs(object);
    nParamPerSeries <- nparam(object);
    structure(object$logLik,nobs=obs,df=nParamPerSeries,class="logLik");
}

#' @importFrom greybox AICc actuals
#' @export
AICc.legion <- function(object, ...){
    llikelihood <- logLik(object);
    llikelihood <- llikelihood[1:length(llikelihood)];
    nSeries <- ncol(actuals(object));
    # Remove covariances in the number of parameters
    nParamAll <- nparam(object) / nSeries - switch(object$loss,
                                                   "likelihood" = nSeries*(nSeries+1)/2,
                                                   "trace" = ,
                                                   "diagonal" = 1);

    obs <- nobs(object);
    IC <- -2*llikelihood + ((2*obs*(nParamAll*nSeries + nSeries*(nSeries+1)/2)) /
                                (obs - (nParamAll + nSeries + 1)));

    return(IC);
}

#' @importFrom greybox BICc
#' @export
BICc.legion <- function(object, ...){
    llikelihood <- logLik(object);
    llikelihood <- llikelihood[1:length(llikelihood)];
    nSeries <- ncol(actuals(object));
    # Remove covariances in the number of parameters
    nParamAll <- nparam(object) / nSeries - switch(object$loss,
                                                   "likelihood" = nSeries*(nSeries+1)/2,
                                                   "trace" = ,
                                                   "diagonal" = 1);

    obs <- nobs(object);
    IC <- -2*llikelihood + (((nParamAll + nSeries*(nSeries+1)/2) *
                                 log(obs * nSeries) * obs * nSeries) /
                                (obs * nSeries - nParamAll - nSeries*(nSeries+1)/2));

    return(IC);
}

#' @importFrom stats nobs
#' @export
nobs.legion <- function(object, ...){
    return(nrow(object$fitted));
}
#' @export
nobs.oves <- function(object, ...){
    return(nrow(object$fitted));
}

#' @importFrom greybox nparam
#' @export
nparam.oves <- function(object, ...){
    nParamReturn <- object$nParam[1,4];
    return(nParamReturn);
}

#' @importFrom stats sigma
#' @export
sigma.legion <- function(object, ...){
    return(object$Sigma);
}

#### Extraction of parameters of models ####
#' @importFrom greybox errorType
#' @export
errorType.legion <- function(object, ...){
    if(any(substr(modelType(object),1,1)==c("A"))){
        return("A");
    }
    else if(any(substr(modelType(object),1,1)==c("M"))){
        return("M");
    }
}

#' @export
coef.legion <- function(object, ...){

    parameters <- object$B;

    return(parameters);
}

#' @importFrom smooth modelType
#' @export
modelType.legion <- function(object, ...){
    model <- object$model;
    modelType <- NA;
    if(!is.null(model)){
        if(gregexpr("VES",model)!=-1){
            modelType <- substring(model,unlist(gregexpr("\\(",model))+1,unlist(gregexpr("\\)",model))-1);
        }
    }

    return(modelType);
}

#### Plotting things ####
#' @export
plot.oves <- function(x, ...){
    ellipsis <- list(...);
    occurrence <- x$occurrence
    if(occurrence=="f"){
        occurrence <- "Fixed probability";
    }
    else if(occurrence=="l"){
        occurrence <- "Logistic probability";
    }
    else{
        occurrence <- "None";
    }

    y <- actuals(x);
    yForecast <- x$forecast;
    yFitted <- x$fitted;
    dataDeltat <- deltat(y);
    forecastStart <- start(yForecast);
    h <- nrow(yForecast);
    nSeries <- ncol(y);
    modelname <- paste0("iVES(",x$model,")")

    pages <- ceiling(nSeries / 5);
    parDefault <- par(no.readonly=TRUE);
    for(j in 1:pages){
        par(mfcol=c(min(5,floor(nSeries/j)),1));
        for(i in 1:nSeries){
            plotRange <- range(min(y[,i],yForecast[,i],yFitted[,i]),
                               max(y[,i],yForecast[,i],yFitted[,i]));
            plot(y[,i],main=paste0(modelname,", series ", i),ylab="Y",
                 ylim=plotRange, xlim=range(time(y[,i])[1],time(yForecast)[max(h,1)]),
                 type="l");
            lines(yFitted[,i],col="purple",lwd=2,lty=2);
            if(h>1){
                lines(yForecast[,i],col="blue",lwd=2);
            }
            else{
                points(yForecast[,i],col="blue",lwd=2,pch=4);
            }
            abline(v=dataDeltat*(forecastStart[2]-2)+forecastStart[1],col="red",lwd=2);
        }
        par(parDefault);
    }
}

#' @export
plot.legion.sim <- function(x, ...){
    ellipsis <- list(...);
    if(is.null(ellipsis$main)){
        ellipsis$main <- x$model;
    }

    if(length(dim(x$data))==2){
        nsim <- 1;
    }
    else{
        nsim <- dim(x$data)[3];
    }

    nSeries <- dim(x$data)[2];
    if(nSeries>10){
        warning("You have generated more than ten time series. We will plot only first ten of them.",
                call.=FALSE);
        x$data <- x$data[,1:10,];
    }

    if(nsim==1){
        if(is.null(ellipsis$ylab)){
            ellipsis$ylab <- "Data";
        }
        ellipsis$x <- x$data;
        do.call(plot, ellipsis);
    }
    else{
        randomNumber <- ceiling(runif(1,1,nsim));
        message(paste0("You have generated ",nsim," time series. Not sure which of them to plot.\n",
                       "Please use plot(ourSimulation$data[,k]) instead. Plotting randomly selected series N",randomNumber,"."));
        if(is.null(ellipsis$ylab)){
            ellipsis$ylab <- paste0("Series N",randomNumber);
        }
        ellipsis$x <- ts(x$data[,,randomNumber]);
        do.call(plot, ellipsis);
    }
}

#### Prints of vector functions ####
#' @export
print.oves <- function(x, ...){

    if(x$probability=="i"){
        occurrence <- "Independent ";
    }
    else if(x$probability=="d"){
        occurrence <- "Dependent ";
    }

    if(x$occurrence=="l"){
        occurrence <- paste0(occurrence,"logistic probability");
    }
    else if(x$occurrence=="f"){
        occurrence <- paste0(occurrence,"fixed probability");
    }
    else{
        occurrence <- "None";
    }
    ICs <- round(c(AIC(x),AICc(x),BIC(x),BICc(x)),4);
    names(ICs) <- c("AIC","AICc","BIC","BICc");
    cat(paste0("Occurrence state space model estimated: ",occurrence,"\n"));
    if(!is.null(x$model)){
        cat(paste0("Underlying ETS model: ",x$model,"\n"));
    }
    cat("Information criteria: \n");
    print(ICs);
}

#' @export
print.legion <- function(x, ...){
    ellipsis <- list(...);
    if(!any(names(ellipsis)=="digits")){
        digits <- 4;
    }
    else{
        digits <- ellipsis$digits;
    }

    holdout <- any(!is.na(x$holdout));
    interval <- any(!is.na(x$PI));

    # if(all(holdout,interval)){
    #     insideinterval <- sum((x$holdout <= x$upper) & (x$holdout >= x$lower)) / length(x$forecast) * 100;
    # }
    # else{
    #     insideinterval <- NULL;
    # }

    intervalType <- x$interval;

    cat(paste0("Time elapsed: ",round(as.numeric(x$timeElapsed,units="secs"),digits)," seconds\n"));
    cat(paste0("Model estimated: ",x$model,"\n"));
    if(!is.null(x$occurrence)){
        if(x$occurrence$probability=="i"){
            occurrence <- "Independent ";
        }
        else if(x$occurrence$probability=="d"){
            occurrence <- "Dependent ";
        }

        if(x$occurrence$occurrence=="l"){
            occurrence <- paste0(occurrence,"logistic probability");
        }
        else if(x$occurrence$occurrence=="f"){
            occurrence <- paste0(occurrence,"fixed probability");
        }
        else{
            occurrence <- "None";
        }

        cat(paste0("Occurrence model estimated: ",occurrence,"\n"));
        if(!is.null(x$occurrence$model)){
            cat(paste0("Underlying ETS model: ",x$model,"\n"));
        }
    }

    cat(paste0("\nLoss function type: ",x$loss))
    if(!is.null(x$lossValue)){
        cat(paste0("; Loss function value: ",round(x$lossValue,digits),"\n"));
    }
    else{
        cat("\n");
    }

    cat("Sample size: "); cat(nobs(x));
    cat("\n");

    if(!is.null(x$nParam)){
        cat("Number of estimated parameters: "); cat(nparam(x));
        cat("\n");

        if(x$nParam[2,4]>0){
            cat("Number of provided parameters: "); cat(x$nParam[2,4]);
            cat("\n");
        }

        cat("Number of series: "); cat(ncol(actuals(x)));
        cat("\n");

        cat("Number of degrees of freedom per series: "); cat(round(nobs(x)-nparam(x) / ncol(actuals(x)),digits));
        cat("\n");
    }

    cat("Information criteria:\n");
    print(round(x$ICs,digits));

    if(interval){
        if(x$interval=="c"){
            intervalType <- "conditional";
        }
        else if(x$interval=="u"){
            intervalType <- "unconditional";
        }
        else if(x$interval=="i"){
            intervalType <- "independent";
        }
        else if(x$interval=="l"){
            intervalType <- "likelihood-based";
        }
        cat(paste0("\n",x$level*100,"% ",intervalType," prediction interval was constructed\n"));
    }

}

#### Simulate data using provided vector object ####
#' @importFrom stats simulate
#' @export
simulate.legion <- function(object, nsim=1, seed=NULL, obs=NULL, ...){
    ellipsis <- list(...);
    if(is.null(obs)){
        obs <- nobs(object);
    }
    if(!is.null(seed)){
        set.seed(seed);
    }

    # Start a list of arguments
    args <- vector("list",0);

    args$nSeries <- ncol(actuals(object));

    if(!is.null(ellipsis$randomizer)){
        randomizer <- ellipsis$randomizer;
    }
    else{
        randomizer <- "rnorm";
    }

    if(randomizer=="rnorm"){
        if(!is.null(ellipsis$mean)){
            args$mean <- ellipsis$mean;
        }
        else{
            args$mean <- 0;
        }

        if(!is.null(ellipsis$sd)){
            args$sd <- ellipsis$sd;
        }
        else{
            args$sd <- sqrt(mean(residuals(object)^2));
        }
    }
    else if(randomizer=="rlaplace"){
        if(!is.null(ellipsis$mu)){
            args$mu <- ellipsis$mu;
        }
        else{
            args$mu <- 0;
        }

        if(!is.null(ellipsis$b)){
            args$b <- ellipsis$b;
        }
        else{
            args$b <- mean(abs(residuals(object)));
        }
    }
    else if(randomizer=="rs"){
        if(!is.null(ellipsis$mu)){
            args$mu <- ellipsis$mu;
        }
        else{
            args$mu <- 0;
        }

        if(!is.null(ellipsis$b)){
            args$b <- ellipsis$b;
        }
        else{
            args$b <- mean(sqrt(abs(residuals(object))));
        }
    }
    else if(randomizer=="mvrnorm"){
        if(!is.null(ellipsis$mu)){
            args$mu <- ellipsis$mu;
        }
        else{
            args$mu <- 0;
        }

        if(!is.null(ellipsis$Sigma)){
            args$Sigma <- ellipsis$Sigma;
        }
        else{
            args$Sigma <- sigma(object);
        }
    }

    args$randomizer <- randomizer;
    args$frequency <- frequency(actuals(object));
    args$obs <- obs;
    args$nsim <- nsim;
    args$initial <- object$initial;

    if(gregexpr("VES",object$model)!=-1){
        if(all(object$phi==1)){
            phi <- 1;
        }
        else{
            phi <- object$phi;
        }
        model <- modelType(object);
        args <- c(args,list(model=model, phi=phi, persistence=object$persistence,
                            transition=object$transition,
                            initialSeason=object$initialSeason));

        simulatedData <- do.call("sim.ves",args);
    }
    else{
        model <- substring(object$model,1,unlist(gregexpr("\\(",object$model))[1]-1);
        message(paste0("Sorry, but simulate is not yet available for the model ",model,"."));
        simulatedData <- NA;
    }
    return(simulatedData);
}

#### Summary of objects ####
#' @export
summary.legion <- function(object, ...){
    print(object);
}
#' @export
summary.oves <- function(object, ...){
    print(object);
}
