association.test.parallel <- function(basename, Y, X, method = c("lm", "lmm"), response = c("quantitative", "binary"), 
                                      test = c("score", "wald", "lrt"), K, eigenK, p = 0, 
                                      tol = .Machine$double.eps^0.25, n.cores = 1, ...) {
  arg <- as.list(match.call())[-1]
  arg$n.cores <- NULL 
  end <- nb.snps(basename)

  if(n.cores == 1) {
    arg$beg = 1;
    arg$end = end;
    return(do.call(association.test.1, arg))
  }
  if(.Platform$OS.type != "unix") {
    stop("FORK cluster unavailable")
  }
  a <- round(seq(1, end, length = n.cores + 1))
  BEG <- a[1:n.cores]
  END <- a[-1]-1;
  END[n.cores] <- end;

  ARG <- rep( list(arg), n.cores)
  ARG <- mapply(function(a, be, en) {a$beg = be; a$end = en; a}, ARG, BEG, END, SIMPLIFY=FALSE )
  cl <- parallel::makeForkCluster(n.cores)
  xx <- parallel::clusterApply(cl, ARG, function(a) do.call(association.test.1, a))
  parallel::stopCluster(cl)
  Reduce(rbind, xx)
}

association.test.1 <- function(basename, Y, X, method = c("lm", "lmm"), response = c("quantitative", "binary"), 
                               test = c("score", "wald", "lrt"), K, eigenK, beg, end, p = 0, 
                               tol = .Machine$double.eps^0.25, ...) {

  if(missing(end))
    x <- read.bed.matrix.part(basename, beg = beg)
  else 
    x <- read.bed.matrix.part(basename, beg = beg, end = end)

  
  arg <- as.list(match.call())[-1]
  arg$beg <- arg$end <- NULL
  arg$basename <- NULL
  arg$x <- x
  do.call(association.test, arg)
}
