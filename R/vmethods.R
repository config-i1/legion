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
    # Number of parameters per series
    nParamPerSeries <- nparam(object)/nSeries - switch(object$loss,
                                                       "likelihood" = (nSeries+1)/2,
                                                       "trace" = ,
                                                       "diagonal" = 1);
    # All the estimated parameters (would differ depending on loss)
    nParamAll <- nparam(object);

    obs <- nobs(object);
    if(obs - (nParamPerSeries + nSeries + 1) <=0){
        IC <- Inf;
    }
    else{
        IC <- -2*llikelihood + 2*(obs*nParamAll / (obs - (nParamPerSeries + nSeries + 1)));
    }

    return(IC);
}

#' @importFrom greybox BICc
#' @export
BICc.legion <- function(object, ...){
    llikelihood <- logLik(object);
    llikelihood <- llikelihood[1:length(llikelihood)];
    nSeries <- ncol(actuals(object));
    # Number of parameters per series
    nParamPerSeries <- nparam(object)/nSeries - switch(object$loss,
                                                       "likelihood" = (nSeries+1)/2,
                                                       "trace" = ,
                                                       "diagonal" = 1);
    # All the estimated parameters (would differ depending on loss)
    nParamAll <- nparam(object);

    obs <- nobs(object);
    if(obs - (nParamPerSeries + nSeries + 1) <=0){
        IC <- Inf;
    }
    else{
        IC <- -2*llikelihood + log(obs)*(obs*nParamAll / (obs - (nParamPerSeries + nSeries + 1)));
    }

    return(IC);
}

#' @importFrom greybox actuals
#' @export
actuals.legion <- function(object, ...){
    return(object$y);
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
    return(object$nParam[1,4]);
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
    else{
        return(NA);
    }
}

#' @export
coef.legion <- function(object, ...){
    return(object$B);
}

#' @importFrom smooth modelType
#' @export
modelType.legion <- function(object, ...){
    model <- object$model;
    modelType <- NA;
    if(!is.null(model)){
        if(gregexpr("VES",model)!=-1 || gregexpr("VETS",model)!=-1){
            modelType <- substring(model,unlist(gregexpr("\\(",model))+1,unlist(gregexpr("\\)",model))-1)[1];
        }
    }

    return(modelType);
}

#### Plotting things ####
#' Plots for the fit and states
#'
#' The function produces diagnostics plots for a \code{legion} model
#'
#' The list of produced plots includes:
#' \enumerate{
#' \item Actuals vs Fitted values. Allows analysing, whether there are any issues in the fit.
#' Does the variability of actuals increase with the increase of fitted values? Is the relation
#' well captured? They grey line on the plot corresponds to the perfect fit of the model.
#' \item Standardised residuals vs Fitted. Plots the points and the confidence bounds
#' (red lines) for the specified confidence \code{level}. Useful for the analysis of outliers;
#' \item Studentised residuals vs Fitted. This is similar to the previous plot, but with the
#' residuals divided by the scales with the leave-one-out approach. Should be more sensitive
#' to outliers;
#' \item Absolute residuals vs Fitted. Useful for the analysis of heteroscedasticity;
#' \item Squared residuals vs Fitted - similar to (3), but with squared values;
#' \item Q-Q plot with the specified distribution. Can be used in order to see if the
#' residuals follow the assumed distribution. The type of distribution depends on the one used
#' in the estimation (see \code{distribution} parameter in \link[greybox]{alm});
#' \item ACF of the residuals. Are the residuals autocorrelated? See \link[stats]{acf} for
#' details;
#' \item Fitted over time. Plots actuals (black line), fitted values (purple line), point forecast
#' (blue line) and prediction interval (grey lines). Can be used in order to make sure that the model
#' did not miss any important events over time;
#' \item Standardised residuals vs Time. Useful if you want to see, if there is autocorrelation or
#' if there is heteroscedasticity in time. This also shows, when the outliers happen;
#' \item Studentised residuals vs Time. Similar to previous, but with studentised residuals;
#' \item PACF of the residuals. No, really, are they autocorrelated? See pacf function from stats
#' package for details;
#' \item Plot of the states of the model. It is not recommended to produce this plot together with
#' the others, because there might be several states, which would cause the creation of a different
#' canvas. In case of "msdecompose", this will produce the decomposition of the series into states
#' on a different canvas.
#' }
#' Which of the plots to produce, is specified via the \code{which} parameter.
#' Currently only \code{which=c(1,4:7)} are supported.
#'
#' @param x Estimated legion model.
#' @param which Which of the plots to produce. The possible options (see details for explanations):
#' \enumerate{
#' \item Actuals vs Fitted values;
#' \item Standardised residuals vs Fitted;
#' \item Studentised residuals vs Fitted;
#' \item Absolute residuals vs Fitted;
#' \item Squared residuals vs Fitted;
#' \item Q-Q plot with the specified distribution;
#' \item Fitted over time;
#' \item Standardised residuals vs Time;
#' \item Studentised residuals vs Time;
#' \item ACF of the residuals;
#' \item PACF of the residuals.
#' \item Plot of states of the model.
#' }
#' @param level Confidence level. Defines width of confidence interval. Used in plots (2), (3), (7), (8),
#' (9), (10) and (11).
#' @param legend If \code{TRUE}, then the legend is produced on plots (2), (3) and (7).
#' @param ask Logical; if \code{TRUE}, the user is asked to press Enter before each plot.
#' @param lowess Logical; if \code{TRUE}, LOWESS lines are drawn on scatterplots, see \link[stats]{lowess}.
#' @param ... The parameters passed to the plot functions. Recommended to use with separate plots.
#' @return The function produces the number of plots, specified in the parameter \code{which}.
#'
#' @template ssAuthor
#' @seealso \link[greybox]{plot.greybox}
#' @keywords ts univar
#' @examples
#'
#' ourModel <- es(c(rnorm(50,100,10),rnorm(50,120,10)), "ANN", h=10)
#' par(mfcol=c(3,4))
#' plot(ourModel, c(1:11))
#' plot(ourModel, 12)
#'
#' @importFrom greybox nvariate is.occurrence graphmaker
#' @importFrom grDevices dev.interactive devAskNewPage
#' @importFrom stats fitted qqline qqnorm
#' @export
plot.legion <- function(x, which=c(1,2,4,6), level=0.95, legend=FALSE,
                        ask=prod(par("mfcol")) < length(which) * nvariate(x) && dev.interactive(),
                        lowess=TRUE, ...){
    nSeries <- nvariate(x);

    # Define, whether to wait for the hit of "Enter"
    if(ask){
        oask <- devAskNewPage(TRUE);
        on.exit(devAskNewPage(oask));
    }

    # 1. Fitted vs Actuals values
    plot1 <- function(x, y, yFitted, i, ...){
        ellipsis <- list(...);

        ellipsis$y <- y;
        ellipsis$x <- yFitted;

        # If this is a mixture model, remove zeroes
        if(is.occurrence(x$occurrence)){
            ellipsis$x <- ellipsis$x[ellipsis$y!=0];
            ellipsis$y <- ellipsis$y[ellipsis$y!=0];
        }

        # Remove NAs
        if(any(is.na(ellipsis$x))){
            ellipsis$y <- ellipsis$y[!is.na(ellipsis$x)];
            ellipsis$x <- ellipsis$x[!is.na(ellipsis$x)];
        }
        if(any(is.na(ellipsis$y))){
            ellipsis$x <- ellipsis$x[!is.na(ellipsis$y)];
            ellipsis$y <- ellipsis$y[!is.na(ellipsis$y)];
        }

        # Title
        if(!any(names(ellipsis)=="main")){
            ellipsis$main <- paste0("Series ",i,", Actuals vs Fitted");
        }
        # If type and ylab are not provided, set them...
        if(!any(names(ellipsis)=="type")){
            ellipsis$type <- "p";
        }
        if(!any(names(ellipsis)=="ylab")){
            ellipsis$ylab <- "Actuals";
        }
        if(!any(names(ellipsis)=="xlab")){
            ellipsis$xlab <- "Fitted";
        }
        # xlim and ylim
        if(!any(names(ellipsis)=="xlim")){
            ellipsis$xlim <- range(c(ellipsis$x,ellipsis$y));
        }
        if(!any(names(ellipsis)=="ylim")){
            ellipsis$ylim <- range(c(ellipsis$x,ellipsis$y));
        }

        # Start plotting
        do.call(plot,ellipsis);
        abline(a=0,b=1,col="grey",lwd=2,lty=2)
        if(lowess){
            lines(lowess(ellipsis$x, ellipsis$y), col="red");
        }
    }

    # 4 and 5. Fitted vs |Residuals| or Fitted vs Residuals^2
    plot3 <- function(x, yResid, yFitted, i,  type="abs", ...){
        ellipsis <- list(...);

        ellipsis$x <- yFitted;
        if(type=="abs"){
            ellipsis$y <- abs(yResid);
        }
        else{
            ellipsis$y <- yResid^2;
        }

        if(is.occurrence(x$occurrence)){
            ellipsis$x <- ellipsis$x[ellipsis$y!=0];
            ellipsis$y <- ellipsis$y[ellipsis$y!=0];
        }
        # Remove NAs
        if(any(is.na(ellipsis$x))){
            ellipsis$x <- ellipsis$x[!is.na(ellipsis$x)];
            ellipsis$y <- ellipsis$y[!is.na(ellipsis$y)];
        }

        if(!any(names(ellipsis)=="main")){
            if(type=="abs"){
                ellipsis$main <- paste0("Series ",i,", |Residuals| vs Fitted");
            }
            else{
                ellipsis$main <- paste0("Series ",i,", Residuals^2 vs Fitted");
            }
        }

        if(!any(names(ellipsis)=="xlab")){
            ellipsis$xlab <- "Fitted";
        }
        if(!any(names(ellipsis)=="ylab")){
            if(type=="abs"){
                ellipsis$ylab <- "|Residuals|";
            }
            else{
                ellipsis$ylab <- "Residuals^2";
            }
        }

        do.call(plot,ellipsis);
        abline(h=0, col="grey", lty=2);
        if(lowess){
            lines(lowess(ellipsis$x, ellipsis$y), col="red");
        }
    }

    # 6. Q-Q with the specified distribution
    plot4 <- function(x, yResid, i, ...){
        ellipsis <- list(...);

        ellipsis$y <- yResid;
        if(is.occurrence(x$occurrence)){
            ellipsis$y <- ellipsis$y[actuals(x$occurrence)!=0];
        }

        if(!any(names(ellipsis)=="xlab")){
            ellipsis$xlab <- "Theoretical Quantile";
        }
        if(!any(names(ellipsis)=="ylab")){
            ellipsis$ylab <- "Actual Quantile";
        }

        if(!any(names(ellipsis)=="main")){
            ellipsis$main <- paste0("Series ",i,", QQ plot of normal distribution");
        }

        do.call(qqnorm, ellipsis);
        qqline(ellipsis$y);
    }

    # 7. Basic plot over time
    plot5 <- function(x, y, yFitted, yHoldout, yForecast, yLower=NULL, yUpper=NULL, ...){
        ellipsis <- list(...);

        ellipsis$actuals <- y;
        if(!is.null(yHoldout)){
            if(is.zoo(ellipsis$actuals)){
                ellipsis$actuals <- zoo(c(as.vector(y),as.vector(yHoldout)),
                                        order.by=c(time(y),time(yHoldout)));
            }
            else{
                ellipsis$actuals <- ts(c(y,yHoldout),
                                       start=start(y),
                                       frequency=frequency(y));
            }
        }
        if(is.null(ellipsis$main)){
            ellipsis$main <- paste0("Series ",i,", ",x$model);
        }
        ellipsis$forecast <- yForecast;
        ellipsis$fitted <- yFitted;
        ellipsis$legend <- FALSE;
        ellipsis$parReset <- FALSE;
        if(!any(x$interval==c("none","n"))){
            ellipsis$lower <- yLower;
            ellipsis$upper <- yUpper;
            ellipsis$level <- x$level;
        }

        do.call(graphmaker, ellipsis);
    }

    # Do plots
    if(any(which==1)){
        for(i in 1:nSeries){
            plot1(x, as.vector(actuals(x)[,i]), as.vector(fitted(x)[,i]), i=i, ...);
        }
    }

    if(any(which %in% c(2,3,8:12))){
        warning("The plots 2, 3, 8-12 are not available yet for the legion class.",
                call.=FALSE);
    }

    if(any(which==4)){
        for(i in 1:nSeries){
            plot3(x, as.vector(residuals(x)[,i]), as.vector(fitted(x)[,i]), i=i, ...);
        }
    }

    if(any(which==5)){
        for(i in 1:nSeries){
            plot3(x, as.vector(residuals(x)[,i]), as.vector(fitted(x)[,i]), i=i, type="squared", ...);
        }
    }

    if(any(which==6)){
        for(i in 1:nSeries){
            plot4(x, as.vector(residuals(x)[,i]), i=i, ...);
        }
    }

    if(any(which==7)){
        for(i in 1:nSeries){
            plot5(x, actuals(x)[,i], fitted(x)[,i], x$holdout, x$forecast[,i], x$lower[,i], x$upper[,i], ...);
        }
    }
}

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
    yFittedted <- x$fitted;
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
            plotRange <- range(min(y[,i],yForecast[,i],yFittedted[,i]),
                               max(y[,i],yForecast[,i],yFittedted[,i]));
            plot(y[,i],main=paste0(modelname,", series ", i),ylab="Y",
                 ylim=plotRange, xlim=range(time(y[,i])[1],time(yForecast)[max(h,1)]),
                 type="l");
            lines(yFittedted[,i],col="purple",lwd=2,lty=2);
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

        cat("Number of series: "); cat(nvariate(x));
        cat("\n");

        cat("Number of degrees of freedom per series: "); cat(round(nobs(x)-nparam(x) / nvariate(x),digits));
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

#### Residuals ####
#' @export
residuals.legion <- function(object, ...){
    return(object$residuals);
}