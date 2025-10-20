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
    return(object$data);
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
    ellipsis <- list(...);
    model <- object$model;
    modelTypeReturned <- NA;
    # If a person asked for PIC part, extract it
    if(!is.null(ellipsis$pic) && ellipsis$pic){
        if(!is.null(model) && gregexpr("VETS",model)!=-1){
            modelTypeReturned <- substring(model,unlist(gregexpr("\\(",model))+1,unlist(gregexpr("\\)",model))-1)[2];
        }
        # Return just basic type for VES
        else{
            modelTypeReturned <- modelType(object);
        }
    }
    else{
        if(!is.null(model)){
            if(gregexpr("VES",model)!=-1 || gregexpr("VETS",model)!=-1){
                modelTypeReturned <- substring(model,unlist(gregexpr("\\(",model))+1,unlist(gregexpr("\\)",model))-1)[1];
            }
        }
    }

    return(modelTypeReturned);
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
#' plot(ourModel, c(1:11))
#' plot(ourModel, 12)
#'
#' @importFrom graphics text
#' @importFrom greybox nvariate is.occurrence graphmaker
#' @importFrom grDevices dev.interactive devAskNewPage
#' @importFrom stats fitted qqline qqnorm
#' @importFrom stats na.pass sd acf pacf qnorm
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

    # 2 and 3: Standardised  / studentised residuals vs Fitted
    plot2 <- function(x, yFitted, yResid, yName, statistic, i, ...){
        ellipsis <- list(...);

        ellipsis$x <- yFitted;
        ellipsis$y <- yResid;

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

        # Main, labs etc
        if(!any(names(ellipsis)=="main")){
            if(errorType(x)=="M"){
                ellipsis$main <- paste0("Series ",i,", log(",yName," Residuals) vs Fitted");
            }
            else{
                ellipsis$main <- paste0("Series ",i,", ",yName," Residuals vs Fitted");
            }
        }

        if(!any(names(ellipsis)=="xlab")){
            ellipsis$xlab <- "Fitted";
        }
        if(!any(names(ellipsis)=="ylab")){
            ellipsis$ylab <- paste0(yName," Residuals");
        }

        if(legend){
            if(ellipsis$x[length(ellipsis$x)]>mean(ellipsis$x)){
                legendPosition <- "bottomright";
            }
            else{
                legendPosition <- "topright";
            }
        }

        outliers <- which(ellipsis$y >statistic[2] | ellipsis$y <statistic[1]);
        # cat(paste0(round(length(outliers)/length(ellipsis$y),3)*100,"% of values are outside the bounds\n"));

        if(!any(names(ellipsis)=="ylim")){
            ellipsis$ylim <- range(c(ellipsis$y,statistic), na.rm=TRUE)*1.2;
            if(legend){
                if(legendPosition=="bottomright"){
                    ellipsis$ylim[1] <- ellipsis$ylim[1] - 0.2*diff(ellipsis$ylim);
                }
                else{
                    ellipsis$ylim[2] <- ellipsis$ylim[2] + 0.2*diff(ellipsis$ylim);
                }
            }
        }

        xRange <- range(ellipsis$x, na.rm=TRUE);
        xRange[1] <- xRange[1] - sd(ellipsis$x, na.rm=TRUE);
        xRange[2] <- xRange[2] + sd(ellipsis$x, na.rm=TRUE);

        do.call(plot,ellipsis);
        abline(h=0, col="grey", lty=2);
        polygon(c(xRange,rev(xRange)),c(statistic[1],statistic[1],statistic[2],statistic[2]),
                col="lightgrey", border=NA, density=10);
        abline(h=statistic, col="red", lty=2);
        if(length(outliers)>0){
            points(ellipsis$x[outliers], ellipsis$y[outliers], pch=16);
            text(ellipsis$x[outliers], ellipsis$y[outliers], labels=outliers, pos=(ellipsis$y[outliers]>0)*2+1);
        }
        if(lowess){
            lines(lowess(ellipsis$x[!is.na(ellipsis$y)], ellipsis$y[!is.na(ellipsis$y)]), col="red");
        }

        if(legend){
            if(lowess){
                legend(legendPosition,
                       legend=c(paste0(round(level,3)*100,"% bounds"),"outside the bounds","LOWESS line"),
                       col=c("red", "black","red"), lwd=c(1,NA,1), lty=c(2,1,1), pch=c(NA,16,NA));
            }
            else{
                legend(legendPosition,
                       legend=c(paste0(round(level,3)*100,"% bounds"),"outside the bounds"),
                       col=c("red", "black"), lwd=c(1,NA), lty=c(2,1), pch=c(NA,16));
            }
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

    # 8 and 9. Standardised / Studentised residuals vs time
    plot6 <- function(x, yResid, yName, statistic, i, ...){

        ellipsis <- list(...);

        ellipsis$x <- yResid;

        if(is.occurrence(x$occurrence)){
            ellipsis$x[ellipsis$x==0] <- NA;
        }

        if(!any(names(ellipsis)=="main")){
            ellipsis$main <- paste0("Series ",i,", ",yName," Residuals vs Time");
        }

        if(!any(names(ellipsis)=="xlab")){
            ellipsis$xlab <- "Time";
        }
        if(!any(names(ellipsis)=="ylab")){
            ellipsis$ylab <- paste0(yName," Residuals");
        }

        # If type and ylab are not provided, set them...
        if(!any(names(ellipsis)=="type")){
            ellipsis$type <- "l";
        }

        outliers <- which(ellipsis$x >statistic[2] | ellipsis$x <statistic[1]);

        if(!any(names(ellipsis)=="ylim")){
            ellipsis$ylim <- c(-max(abs(ellipsis$x),na.rm=TRUE),
                               max(abs(ellipsis$x),na.rm=TRUE))*1.2;
        }

        if(legend){
            legendPosition <- "topright";
            ellipsis$ylim[2] <- ellipsis$ylim[2] + 0.2*diff(ellipsis$ylim);
            ellipsis$ylim[1] <- ellipsis$ylim[1] - 0.2*diff(ellipsis$ylim);
        }

        # Start plotting
        do.call(plot,ellipsis);
        if(length(outliers)>0){
            points(time(ellipsis$x)[outliers], ellipsis$x[outliers], pch=16);
            text(time(ellipsis$x)[outliers], ellipsis$x[outliers], labels=outliers, pos=(ellipsis$x[outliers]>0)*2+1);
        }
        # If there is occurrence model, plot points to fill in breaks
        if(is.occurrence(x$occurrence)){
            points(time(ellipsis$x), ellipsis$x);
        }
        if(lowess){
            # Substitute NAs with the mean
            if(any(is.na(ellipsis$x))){
                ellipsis$x[is.na(ellipsis$x)] <- mean(ellipsis$x, na.rm=TRUE);
            }
            lines(lowess(c(1:length(ellipsis$x)),ellipsis$x), col="red");
        }
        abline(h=0, col="grey", lty=2);
        abline(h=statistic[1], col="red", lty=2);
        abline(h=statistic[2], col="red", lty=2);
        polygon(c(1:nobs(x), c(nobs(x):1)),
                c(rep(statistic[1],nobs(x)), rep(statistic[2],nobs(x))),
                col="lightgrey", border=NA, density=10);
        if(legend){
            legend(legendPosition,legend=c("Residuals",paste0(level*100,"% prediction interval")),
                   col=c("black","red"), lwd=rep(1,3), lty=c(1,1,2));
        }
    }

    # 10 and 11. ACF and PACF
    plot7 <- function(x, yResid, type="acf", i, ...){
        ellipsis <- list(...);

        if(!any(names(ellipsis)=="main")){
            if(type=="acf"){
                ellipsis$main <- paste0("Series ",i,", Autocorrelation Function of Residuals");
            }
            else{
                ellipsis$main <- paste0("Series ",i,", Partial Autocorrelation Function of Residuals");
            }
        }

        if(!any(names(ellipsis)=="xlab")){
            ellipsis$xlab <- "Lags";
        }
        if(!any(names(ellipsis)=="ylab")){
            if(type=="acf"){
                ellipsis$ylab <- "ACF";
            }
            else{
                ellipsis$ylab <- "PACF";
            }
        }

        if(!any(names(ellipsis)=="ylim")){
            ellipsis$ylim <- c(-1,1);
        }

        if(type=="acf"){
            theValues <- acf(yResid, plot=FALSE, na.action=na.pass);
        }
        else{
            theValues <- pacf(yResid, plot=FALSE, na.action=na.pass);
        }
        ellipsis$x <- theValues$acf[-1];
        statistic <- qnorm(c((1-level)/2, (1+level)/2),0,sqrt(1/nobs(x)));

        ellipsis$type <- "h"

        do.call(plot,ellipsis);
        abline(h=0, col="black", lty=1);
        abline(h=statistic, col="red", lty=2);
        if(any(ellipsis$x>statistic[2] | ellipsis$x<statistic[1])){
            outliers <- which(ellipsis$x >statistic[2] | ellipsis$x <statistic[1]);
            points(outliers, ellipsis$x[outliers], pch=16);
            text(outliers, ellipsis$x[outliers], labels=outliers, pos=(ellipsis$x[outliers]>0)*2+1);
        }
    }

    # 12. Plot of states
    plot8 <- function(x, ...){
        ellipsis <- list(...);

        parDefault <- par(no.readonly = TRUE);
        on.exit(par(parDefault));
        statesNames <- c(colnames(x$states),
                         paste0("residuals of ",colnames(actuals(x))));
        x$states <- cbind(x$states,residuals(x));
        colnames(x$states) <- statesNames;
        if(ncol(x$states)>10){
            message("Too many states. Plotting them one by one on several graphs.");
            if(is.null(ellipsis$main)){
                ellipsisMain <- NULL;
            }
            else{
                ellipsisMain <- ellipsis$main;
            }
            nPlots <- ceiling(ncol(x$states)/10);
            for(i in 1:nPlots){
                if(is.null(ellipsisMain)){
                    ellipsis$main <- paste0("States of ",x$model,", part ",i);
                }
                ellipsis$x <- x$states[,(1+(i-1)*10):min(i*10,ncol(x$states)),drop=FALSE];
                do.call(plot, ellipsis);
            }
        }
        else{
            if(ncol(x$states)<=5){
                ellipsis$nc <- 1;
            }
            if(is.null(ellipsis$main)){
                ellipsis$main <- paste0("States of ",x$model);
            }
            ellipsis$x <- x$states;
            do.call(plot, ellipsis);
        }
    }

    # Do plots
    for(j in which){
        if(any(j==1)){
            for(i in 1:nSeries){
                plot1(x, as.vector(actuals(x)[,i]), as.vector(fitted(x)[,i]), i, ...);
            }
        }
        else if(any(j==2)){
            # Get the residuals and statistic for outliers
            errors <- rstandard(x);
            outliers <- outlierdummy(x, level=level, type="rstandard");
            statistic <- outliers$statistic;

            for(i in 1:nSeries){
                plot2(x, as.vector(fitted(x)[,i]), as.vector(errors[,i]), "Standardised", statistic, i, ...);
            }
        }
        else if(any(j==3)){
            # Get the residuals and statistic for outliers
            errors <- rstudent(x);
            outliers <- outlierdummy(x, level=level, type="rstudent");
            statistic <- outliers$statistic;

            for(i in 1:nSeries){
                plot2(x, as.vector(fitted(x)[,i]), as.vector(errors[,i]), "Studentised", statistic, i, ...);
            }
        }
        else if(any(j==4)){
            for(i in 1:nSeries){
                plot3(x, as.vector(residuals(x)[,i]), as.vector(fitted(x)[,i]), i, ...);
            }
        }
        else if(any(j==5)){
            for(i in 1:nSeries){
                plot3(x, as.vector(residuals(x)[,i]), as.vector(fitted(x)[,i]), i, type="squared", ...);
            }
        }
        else if(any(j==6)){
            for(i in 1:nSeries){
                plot4(x, as.vector(residuals(x)[,i]), i, ...);
            }
        }
        else if(any(j==7)){
            for(i in 1:nSeries){
                plot5(x, actuals(x)[,i], fitted(x)[,i], x$holdout[,i], x$forecast[,i], x$PI[,i*2-1], x$PI[,i*2], ...);
            }
        }
        else if(any(j==8)){
            # Get the residuals and statistic for outliers
            errors <- rstandard(x);
            outliers <- outlierdummy(x, level=level, type="rstandard");
            statistic <- outliers$statistic;

            for(i in 1:nSeries){
                plot6(x, errors[,i], "Standardised", statistic, i, ...);
            }
        }
        else if(any(j==9)){
            # Get the residuals and statistic for outliers
            errors <- rstudent(x);
            outliers <- outlierdummy(x, level=level, type="rstudent");
            statistic <- outliers$statistic;

            for(i in 1:nSeries){
                plot6(x, errors[,i], "Studentised", statistic, i, ...);
            }
        }
        else if(any(j==10)){
            errors <- residuals(x);
            for(i in 1:nSeries){
                plot7(x, errors[,i], "acf", i, ...);
            }
        }
        else if(any(j==11)){
            errors <- residuals(x);
            for(i in 1:nSeries){
                plot7(x, errors[,i], "pacf", i, ...);
            }
        }
        else if(any(j==12)){
            plot8(x, ...);
        }
    }
}

#' @export
plot.oves <- function(x, ...){
    ellipsis <- list(...);
    occurrence <- x$occurrence
    if(occurrence=="fixed"){
        occurrence <- "Fixed probability";
    }
    else if(occurrence=="logistic"){
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
    on.exit(par(parDefault));
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

    if(x$probability=="independent"){
        occurrence <- "Independent ";
    }
    else if(x$probability=="dependent"){
        occurrence <- "Dependent ";
    }

    if(x$occurrence=="logistic"){
        occurrence <- paste0(occurrence,"logistic probability");
    }
    else if(x$occurrence=="fixed"){
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
    # interval <- any(!is.na(x$PI));
    #
    # if(all(holdout,interval)){
    #     insideinterval <- sum((x$holdout <= x$upper) & (x$holdout >= x$lower)) / length(x$forecast) * 100;
    # }
    # else{
    #     insideinterval <- NULL;
    # }

    # intervalType <- x$interval;

    cat(paste0("Time elapsed: ",round(as.numeric(x$timeElapsed,units="secs"),digits)," seconds\n"));
    cat(paste0("Model estimated: ",x$model,"\n"));
    if(!is.null(x$occurrence)){
        if(x$occurrence$probability=="independent"){
            occurrence <- "Independent ";
        }
        else if(x$occurrence$probability=="dependent"){
            occurrence <- "Dependent ";
        }

        if(x$occurrence$occurrence=="logistic"){
            occurrence <- paste0(occurrence,"logistic probability");
        }
        else if(x$occurrence$occurrence=="fixed"){
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

    if(x$loss!="custom"){
        cat("Information criteria:\n");
        print(round(x$ICs,digits));
    }
    else{
        cat("Information criteria are not available for the defined loss.\n");
    }

    # if(interval){
    #     if(x$interval=="c"){
    #         intervalType <- "conditional";
    #     }
    #     else if(x$interval=="u"){
    #         intervalType <- "unconditional";
    #     }
    #     else if(x$interval=="i"){
    #         intervalType <- "independent";
    #     }
    #     else if(x$interval=="l"){
    #         intervalType <- "likelihood-based";
    #     }
    #     cat(paste0("\n",x$level*100,"% ",intervalType," prediction interval was constructed\n"));
    # }

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

    args$nvariate <- nvariate(object);

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

# Function decomposes the sigma matrix and returns the Sigma^{-1/2}
# This is an internal function, needed for rstandard and rstudent
sigmasqrtCalculator <- function(sigmaMatrix,nSeries){
    sigmaMatrixEigen <- eigen(sigmaMatrix,symmetric=TRUE);
    return(sigmaMatrixEigen$vectors %*% diag(sigmaMatrixEigen$values^{-0.5}) %*%
               solve(sigmaMatrixEigen$vectors,diag(nSeries)));

}

#' @importFrom stats rstandard
#' @export
rstandard.legion <- function(model, ...){
    obs <- nobs(model);
    nSeries <- nvariate(model);
    rstandardValues <- residuals(model);
    errors <- t(residuals(model));
    # If this is an occurrence model, then only modify the non-zero obs
    # Also, if there are NAs in actuals, consider them as occurrence
    # if(is.occurrence(model$occurrence)){
    #     residsToGo <- which(actuals(model$occurrence)!=0 & !is.na(actuals(model)));
    # }
    # else{
    residsToGo <- c(1:obs);
    # }

    sigmaMatrix <- sigma(model);
    # Decompose the sigma matrix and get the Sigma^{-1/2}
    sigmaMatrixSQRT <- sigmasqrtCalculator(sigmaMatrix,nSeries);

    # Either additive or multiplicative error -> Normal or Log Normal
    if(errorType(model)=="A"){
        rstandardValues[] <- t(sigmaMatrixSQRT %*% (errors - rowMeans(errors[,residsToGo])));
    }
    else{
        # Debias the residuals
        errors[] <- errors - diag(sigmaMatrix);
        rstandardValues[] <- t(sigmaMatrixSQRT %*% errors);
    }
    return(rstandardValues);
}

#' @importFrom stats rstudent
#' @export
rstudent.legion <- function(model, ...){
    obs <- nobs(model);
    nSeries <- nvariate(model);
    df <- obs - nparam(model) / nSeries;
    rstudentised <- errors <- t(residuals(model));
    rstudentValues <- errorsT <- residuals(model);
    # If this is an occurrence model, then only modify the non-zero obs
    # Also, if there are NAs in actuals, consider them as occurrence
    # if(is.occurrence(model$occurrence)){
    #     residsToGo <- which(actuals(model$occurrence)!=0 & !is.na(actuals(model)));
    # }
    # else{
    residsToGo <- c(1:obs);
    # }

    sigmaMatrix <- sigma(model);

    # Either additive or multiplicative error -> Normal or Log Normal
    if(errorType(model)=="A"){
        errors[] <- errors - rowMeans(errors[,residsToGo]);
    }
    else{
        # Debias the residuals
        errors[] <- errors - diag(sigmaMatrix);
    }

    # Go through each observation, remove it and recalculate sigma matrix
    for(i in residsToGo){
        sigmaMatrix[] <- (errors[,-i,drop=FALSE] %*% errorsT[-i,,drop=FALSE]) / df;
        sigmaMatrixSQRT <- sigmasqrtCalculator(sigmaMatrix,nSeries);
        rstudentised[,i] <- sigmaMatrixSQRT %*% errors[,i,drop=FALSE];
    }

    # Return the original object
    rstudentValues[] <- t(rstudentised);

    return(rstudentValues);
}

#' @importFrom greybox outlierdummy
#' @export
outlierdummy.legion <- function(object, level=0.999, type=c("rstandard","rstudent"), ...){
    # Function returns the matrix of dummies with outliers
    type <- match.arg(type);
    errors <- switch(type,"rstandard"=rstandard(object),"rstudent"=rstudent(object));
    statistic <- qnorm(c((1-level)/2, (1+level)/2), 0, 1);
    outliersID <- which(errors>statistic[2] | errors<statistic[1], arr.ind=TRUE);
    outliersNumber <- nrow(outliersID);
    if(outliersNumber>0){
        outliers <- matrix(0, nobs(object), nvariate(object),
                           dimnames=list(rownames(actuals(object)), colnames(actuals(object))));
        outliers[outliersID] <- 1;
    }
    else{
        outliers <- NULL;
    }
    outliersID <- which(errors>statistic[2] | errors<statistic[1], arr.ind=FALSE);

    return(structure(list(outliers=outliers, statistic=statistic, id=outliersID,
                          level=level, type=type),
                     class="outlierdummy"));
}

#### Forecasts ####
#' @importFrom generics forecast
#' @importFrom stats qchisq
#' @export
forecast.legion <- function(object, h=10, #newdata=NULL, occurrence=NULL,
                            interval=c("none", "prediction"),
                            level=0.95, side=c("both","upper","lower"), cumulative=FALSE, nsim=10000, ...){
    # Forecast function for VES and VETS
    interval <- match.arg(interval);
    side <- match.arg(side);

    obsInSample <- nobs(object);
    nSeries <- nvariate(object);
    nParam <- nparam(object);
    lagsModel <- object$lagsAll;
    lagsModelMax <- max(lagsModel);

    matF <- object$transition;
    matG <- object$persistence;
    matW <- object$measurement;
    matVt <- t(object$states[obsInSample+(1:lagsModelMax),,drop=FALSE]);

    model <- modelType(object);
    # If chosen model is "AAdN" or anything like that, we are taking the appropriate values
    Etype <- substring(model,1,1);
    Ttype <- substring(model,2,2);
    Stype <- substring(model,nchar(model),nchar(model));

    yIndex <- time(actuals(object));
    yClasses <- class(actuals(object));
    # Create indices for the future
    if(any(yClasses=="ts")){
        # ts structure
        yFrequency <- frequency(actuals(object));
        if(is.null(object$holdout)){
            yForecastStart <- time(actuals(object))[obsInSample]+deltat(actuals(object));
        }
        else{
            yForecastStart <- time(object$holdout)[1];
        }
    }
    else{
        # zoo
        yIndex <- time(actuals(object));
        if(is.null(object$holdout)){
            yForecastIndex <- yIndex[obsInSample]+diff(tail(yIndex,2))*c(1:h);
        }
        else{
            if(nrow(object$holdout)>=h){
                yForecastIndex <- time(object$holdout)[1:h];
            }
            else{
                yForecastIndex <- c(time(object$holdout), tail(time(object$holdout),1)+
                                    as.numeric(diff(tail(yIndex,2)))*c(1:(nrow(object$holdout)-h)));
            }
        }
    }
    yNames <- colnames(actuals(object));
    if(is.null(yNames)){
        yNames <- paste0("Series",1:nSeries);
    }

    #### Point forecasts ####
    if(any(yClasses=="ts")){
        yForecast <- ts(matrix(NA,h,nSeries,
                               dimnames=list(NULL,yNames)),
                        start=yForecastStart, frequency=yFrequency);
    }
    else{
        yForecast <- zoo(matrix(NA,h,nSeries,
                                dimnames=list(NULL,yNames)),
                         order.by=yForecastIndex);
    }

    yForecast[] <- t(vForecasterWrap(matVt, Matrix(matF, sparse=TRUE), Matrix(matW, sparse=TRUE),
                                     nSeries, h, Etype, Ttype, Stype, lagsModel));

    if(cumulative){
        yForecast <- colSums(yForecast);
    }

    #### Intervals ####
    # Deal with prediction intervals
    Sigma <- sigma(object);
    PI <- NA;
    if(interval!="none"){
        nElements <- length(lagsModel);
        # Number of degrees of freedom per series
        df <- obsInSample - nParam/nSeries;

        # In case of individual we use either Z distribution or Chebyshev inequality
        if(interval=="prediction"){
            if(df>0){
                quantUpper <- qnorm((1+level)/2,0,1);
                quantLower <- qnorm((1-level)/2,0,1);
            }
            else{
                quantUpper <- sqrt(1/((1-level)/2));
                quantLower <- -quantUpper;
            }
        }
        # In case of conditional / unconditional, we use Chi-squared distribution
        else{
            quant <- qchisq(level,df=nSeries);
        }

        nPoints <- 100;
        if(interval=="conditional"){
            # Number of points in the ellipse
            PI <- array(NA, c(h,2*nPoints^(nSeries-1),nSeries),
                        dimnames=list(paste0("h",c(1:h)), NULL,
                                      yNames));
        }
        else{
            PI <- matrix(NA, nrow=h, ncol=nSeries*2,
                         dimnames=list(paste0("h",c(1:h)),
                                       paste0(rep(yNames,each=2),c("_lower","_upper"))));
        }

        # Array of final variance matrices
        varVec <- array(NA,c(h,nSeries,nSeries));
        # This is needed for the first observations, where we do not care about the transition equation
        for(i in 1:min(h,lagsModelMax)){
            varVec[i,,] <- Sigma;
        }

        if(h>1){
            if(cumulative){
                covarVec <- array(NA,c(h,nSeries,nSeries));
            }

            matrixOfVarianceOfStates <- array(0,c(nElements,nElements,h+lagsModelMax));
            # This multiplication does not make sense
            matrixOfVarianceOfStates[,,1:lagsModelMax] <- matG %*% Sigma %*% t(matG);
            matrixOfVarianceOfStatesLagged <- as.matrix(matrixOfVarianceOfStates[,,1]);

            # New transition and measurement for the internal use
            matFNew <- matrix(0,nElements,nElements);
            matWNew <- matrix(0,nSeries,nElements);

            # selectionMat is needed for the correct selection of lagged variables in the array
            # elementsNew are needed for the correct fill in of all the previous matrices
            selectionMat <- matFNew;
            elementsNew <- rep(FALSE,nElements);

            # Define chunks, which correspond to the lags with h being the final one
            chuncksOfHorizon <- c(1,unique(lagsModel),h);
            chuncksOfHorizon <- sort(chuncksOfHorizon);
            chuncksOfHorizon <- chuncksOfHorizon[chuncksOfHorizon<=h];
            chuncksOfHorizon <- unique(chuncksOfHorizon);

            # Length of the vector, excluding the h at the end
            chunksLength <- length(chuncksOfHorizon) - 1;

            elementsNew <- lagsModel<=(chuncksOfHorizon[1]);
            matWNew[,elementsNew] <- matW[,elementsNew];

            for(j in 1:chunksLength){
                selectionMat[lagsModel==chuncksOfHorizon[j],] <- chuncksOfHorizon[j];
                selectionMat[,lagsModel==chuncksOfHorizon[j]] <- chuncksOfHorizon[j];

                elementsNew <- lagsModel < (chuncksOfHorizon[j]+1);
                matFNew[,elementsNew] <- matF[,elementsNew];
                matWNew[,elementsNew] <- matW[,elementsNew];

                for(i in (chuncksOfHorizon[j]+1):chuncksOfHorizon[j+1]){
                    selectionMat[lagsModel>chuncksOfHorizon[j],] <- i;
                    selectionMat[,lagsModel>chuncksOfHorizon[j]] <- i;

                    matrixOfVarianceOfStatesLagged[elementsNew,
                                                   elementsNew] <- matrixOfVarianceOfStates[
                                                       cbind(rep(c(1:nElements),
                                                                 each=nElements),
                                                             rep(c(1:nElements),
                                                                 nElements),
                                                             i - c(selectionMat))];

                    matrixOfVarianceOfStates[,,i] <- (matFNew %*%
                                                          matrixOfVarianceOfStatesLagged %*% t(matFNew) +
                                                          matG %*% Sigma %*% t(matG));
                    varVec[i,,] <- matWNew %*% matrixOfVarianceOfStatesLagged %*% t(matWNew) + Sigma;
                    if(cumulative){
                        covarVec[i] <- matWNew %*% matFNew %*% matG;
                    }
                }
            }

            if(cumulative){
                varVec <- apply(varVec,c(2,3),sum) + 2*Sigma %*% apply(covarVec*array(c(0,h:2),
                                                                                      c(h,nSeries,nSeries)),
                                                                       c(2,3),sum);
            }
        }

        #### Produce PI matrix
        # This one doesn't work yet
        if(any(interval==c("conditional","unconditional"))){
            # eigensList contains eigenvalues and eigenvectors of the covariance matrix
            eigensList <- apply(varVec,1,eigen);
            # eigenLimits specify the lowest and highest ellipse points in all dimensions
            eigenLimits <- matrix(NA,nSeries,2);
            # ellipsePoints contains coordinates of the ellipse on the eigenvectors basis
            ellipsePoints <- array(NA, c(h, 2*nPoints^(nSeries-1), nSeries));
            for(i in 1:h){
                eigenLimits[,2] <- sqrt(quant / eigensList[[i]]$value);
                eigenLimits[,1] <- -eigenLimits[,2];
                ellipsePoints[i,,nSeries] <- rep(seq(eigenLimits[nSeries,1],
                                                     eigenLimits[nSeries,2],
                                                     length.out=nPoints),nSeries);
                for(j in (nSeries-1):1){
                    ellipsePoints[i,,nSeries];
                }
            }
        }
        else if(interval=="prediction"){
            variances <- apply(varVec,1,diag);
            for(i in 1:nSeries){
                PI[,2*i-1] <- quantLower * sqrt(variances[i,]);
                PI[,2*i] <- quantUpper * sqrt(variances[i,]);
            }
        }

        for(i in 1:nSeries){
            PI[,i*2-1] <- PI[,i*2-1] + yForecast[,i];
            PI[,i*2] <- PI[,i*2] + yForecast[,i];
        }

        if(any(yClasses=="ts")){
            PI <-  ts(PI, start=yForecastStart, frequency=yFrequency);
        }
        else{
            PI <-  zoo(PI, order.by=yForecastIndex);
        }
    }

    # Take exponent of the forecasts to get back to the original scale
    if(Etype=="M"){
        yForecast[] <- exp(yForecast);
        PI[] <- exp(PI);
    }

    # This needs to be treated separately via forecast.oves
    # if(occurrence!="n"){
    #     if(!occurrenceModelProvided){
    #         ovesModel <- oves(ts(t(ot),frequency=yFrequency),
    #                           occurrence=occurrence, h=h, holdout=FALSE,
    #                           probability="dependent", model=ovesModel);
    #     }
    #     yForecast[] <- yForecast * t(ovesModel$forecast);
    # }

    return(structure(list(mean=yForecast, PI=PI, model=object,
                          level=level, interval=interval, side=side, cumulative=cumulative, h=h),
                     class=c("legion.forecast","smooth.forecast","forecast")))
}

#' @export
print.legion.forecast <- function(x, ...){
    nSeries <- ncol(x$mean);

    if(x$interval!="none"){
        returnedValue <- switch(x$side,
                                "both"=cbind(x$mean,x$PI),
                                "lower"=cbind(x$mean,x$PI[,c(1:nSeries)*2-1]),
                                "upper"=cbind(x$mean,x$PI[,c(1:nSeries)*2]));
        colnames(returnedValue) <- switch(x$side,
                                          "both"=c(paste0(colnames(x$mean), "_Mean"),colnames(x$PI)),
                                          "lower"=c(paste0(colnames(x$mean), "_Mean"),colnames(x$PI)[c(1:nSeries)*2-1]),
                                          "upper"=c(paste0(colnames(x$mean), "_Mean"),colnames(x$PI)[c(1:nSeries)*2]))
    }
    else{
        returnedValue <- x$mean;
    }
    print(returnedValue);
}

#' @export
plot.legion.forecast <- function(x, ...){
    x$forecast <- x$mean;
    x$holdout <- x$model$holdout;
    x$data <- actuals(x);
    x$fitted <- fitted(x);
    x$model <- x$model$model;
    class(x) <- "legion";
    plot.legion(x, which=7, ...)
}