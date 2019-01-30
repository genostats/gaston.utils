association.test.approx.pql <- function(x, Y = x@ped$pheno, X = matrix(1, nrow(x)), 
                             eigenK, beg = 1, end = ncol(x), p = 0, 
                             tol = .Machine$double.eps^0.25) {

  if(beg < 1 | end > ncol(x)) stop("range too wide")
  if(is.null(x@mu) | is.null(x@p)) stop("Need mu and p to be set in x (use set.stats)")
  if(length(Y) != nrow(x)) stop("Dimensions of Y and x mismatch")
  
  X <- as.matrix(X)
  if(nrow(X) != nrow(x)) stop("Dimensions of Y and x mismatch")
  
  # check dimensions before anything
  n <- nrow(x)
  if(n != nrow(eigenK$vectors) | n != ncol(eigenK$vectors) | n != length(eigenK$values)) 
    stop("eigenK and x dimensions don't match")

  # preparation de X 
  # if(p > 0) {
  #    X <- cbind(X, eigenK$vectors[,seq_len(p)])
  #    X <- gaston:::trans.X(X, mean.y = mean(Y))
  #} else {
  #  X <- gaston:::trans.X(X, mean.y = mean(Y))
  #}


      if(missing(eigenK)) 
        stop("eigenK is mandatory")
      if( any(is.na(Y)) ) 
        stop("Can't handle missing data in Y, please recompute eigenK for the individuals with non-missing phenotype")

        X <- cbind(X, 0) # space for the SNP
        t <- .Call("gg_GWAS_approx_pql_bed", PACKAGE = "gaston.utils", x@bed, x@p, Y, X, eigenK$values, eigenK$vectors, beg-1, end-1, tol)
        t$p <- pchisq( (t$beta/t$sd)**2, df = 1, lower.tail=FALSE)

  L <- data.frame(chr = x@snps$chr, pos = x@snps$pos, id  = x@snps$id,  A1 = x@snps$A1, A2 = x@snps$A2, freqA2 = x@p)
  if(beg > 1 | end < ncol(x))  # avoid copy
    L <- L[beg:end,] 

  data.frame( c( L, t) )
}


