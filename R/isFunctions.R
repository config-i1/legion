#' legion classes checkers
#'
#' Functions to check if an object is of the specified class
#'
#' The list of methods includes:
#' \itemize{
#' \item \code{is.legion()} tests if the object was produced by a vector model (e.g.
#' \link[legion]{ves});
#' \item \code{is.oves()} tests if the object was produced by \link[legion]{oves}
#' function;
#' \item \code{is.legion.sim()} tests if the object was produced by the functions
#' \link[legion]{sim.ves};
#' }
#'
#' @param x The object to check.
#' @return \code{TRUE} if this is the specified class and \code{FALSE} otherwise.
#'
#' @template ssAuthor
#' @keywords ts univar
#' @examples
#'
#' ourModel <- ves(cbind(rnorm(100,100,10),rnorm(100,100,10)))
#'
#' is.legion(ourModel)
#'
#' @rdname isFunctions
#' @export
is.legion <- function(x){
    return(inherits(x,"legion"))
}

#' @rdname isFunctions
#' @export
is.oves <- function(x){
    return(inherits(x,"oves"))
}

#' @rdname isFunctions
#' @export
is.legion.sim <- function(x){
    return(inherits(x,"legion.sim"))
}
