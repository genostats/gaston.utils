approx.pql <- function(Y, X = matrix(1, nrow=length(Y)), eigenK, verbose = getOption("gaston.verbose", TRUE), tol = .Machine$double.eps^0.25) {

  if( any(is.na(Y)) ) 
    stop('Missing data in Y.')

  if( length(Y)!=nrow(X) ) 
    stop('Length of Y and the number of rows of X differ.')

  Sigma <- eigenK$values
  if( length(Y)!=length(Sigma) ) 
    stop('Length of Y and number of eigenvalues differ.')

  w <- which(Sigma < 1e-6)
  Sigma[w] <- 1e-6

  U <- eigenK$vectors
 
  return(.Call('gg_logitmm_approxpql', PACKAGE = "gaston.utils", Y, X, Sigma, U, tol, verbose))
  
}


