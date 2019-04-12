association.test.approx.pql <- function(x, Y = x@ped$pheno, X = matrix(1, nrow(x)), K, beg = 1, end = ncol(x), p = 0, 
                             tol = .Machine$double.eps^0.25, dominance = FALSE, ...) {

  if(beg < 1 | end > ncol(x)) stop("range too wide")
  if(is.null(x@mu) | is.null(x@p)) stop("Need mu and p to be set in x (use set.stats)")
  if(length(Y) != nrow(x)) stop("Dimensions of Y and x mismatch")
  
  X <- as.matrix(X)
  if(nrow(X) != nrow(x)) stop("Dimensions of Y and x mismatch")
  
  # check dimensions before anything
  n <- nrow(x)
  if(missing(K)) stop("argument K is mandatory")
 
  if(!is.list(K)) {
    if(n != nrow(K) | n != ncol(K))
      stop("K and x dimensions don't match")
  } else {
    if(any(n != sapply(K, nrow)) | any(n != sapply(K, ncol)))
      stop("K and x dimensions don't match")
  }

  # preparation de X 
  if(p > 0) {
    X <- cbind(X, eigenK$vectors[,seq_len(p)])
    X <- gaston:::trans.X(X, mean.y = mean(Y))
  } else {
    X <- gaston:::trans.X(X, mean.y = mean(Y))
  }

  # on peut gérer les données manquantes dans Y
  if( any(is.na(Y)) ) {
    w <- !is.na(Y)
    X <- as.matrix(X[w,])
    Y <- Y[w]
    if (is.list(K)) K <- lapply(K, function(x) x[w,w]) else K <- K[w,w]
    warning(sum(!w), 'individuals with missing phenotype are ignored.\n')
  } 
   
  # c'est parti
  model <- logistic.mm.aireml(Y, X = X, K, get.P = TRUE, ... )
  omega <- model$BLUP_omega
  if (!is.null(X)) omega <- omega + X%*%model$BLUP_beta
  pi <- 1/(1+exp(-omega))

  t <- .Call("gg_GWAS_approx_pql_bed", PACKAGE = "gaston.utils", x@bed, Y-pi, model$P, x@p, beg-1, end-1)
  t$p <- pchisq( (t$beta/t$sd)**2, df = 1, lower.tail = FALSE)

  # mise en forme
  L <- data.frame(chr = x@snps$chr, pos = x@snps$pos, id  = x@snps$id,  A1 = x@snps$A1, A2 = x@snps$A2, freqA2 = x@p)
  if(beg > 1 | end < ncol(x))  # avoid copy
    L <- L[beg:end,] 

  data.frame( c( L, t) )
}


