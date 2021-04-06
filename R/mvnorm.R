#' @importFrom stats dnorm
# @export
dmvnorm <- function(q, mean=0, Sigma=1, log=FALSE){
    # The function returns PDF of multivariate normal distribution
    # q should contain obs in columns and series in rows
    if(!is.null(ncol(q))){
        obs <- ncol(q);
        nSeries <- nrow(q);
    }
    else{
        return(dnorm(x=q, mean=mean, sd=sqrt(Sigma), log=log));
    }
    # If dims of mean differ from the q, create the matrix
    # if(!all(dim(q)==dim(mean))){
    #     mean <- matrix(mean,nSeries,obs,byrow=TRUE);
    # }
    # Take invert. If it doesn't work, return NAs
    SigmaInv <- try(chol2inv(chol(Sigma)));
    if(inherits(SigmaInv,"try-error")){
        SigmaInv <- try(solve(Sigma, diag(nSeries), tol=1e-20), silent=TRUE);
        if(inherits(SigmaInv,"try-error")){
            return(rep(NA,obs));
        }
    }
    # Defin X and X transposed
    xt <- q - mean;
    x <- t(xt);
    # Calculate PDF
    mvnormReturn <- vector("numeric",obs);
    for(i in 1:obs){
        mvnormReturn[i] <- x[i,,drop=FALSE] %*% SigmaInv %*% xt[,i,drop=FALSE];
    }
    mvnormReturn[] <- exp(-0.5 * mvnormReturn) * (2*pi)^{-nSeries/2} *det(Sigma)^{-0.5};
    if(log){
        mvnormReturn[] <- log(mvnormReturn);
    }
    return(mvnormReturn);
}