association.test.dosage.parallel <- function(filename, Y, X, method = c("lm", "lmm"), response = c("quantitative", "binary"), 
                                      test = c("score", "wald", "lrt"), K, eigenK, p = 0, 
                                      tol = .Machine$double.eps^0.25, n.cores = 1, ...) {
  nb.snps <- R.utils::countLines(filename)
  nb.inds <- nb.inds.gen.file(filename)

  arg <- as.list(match.call())[-1]
  arg$n.cores <- NULL 
  end <- nb.snps

  if(n.cores == 1) {
    arg$beg = 1;
    arg$end = end;
    return(do.call(association.test.dosage, arg))
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
  xx <- parallel::clusterApply(cl, ARG, function(a) do.call(association.test.dosage, a))
  parallel::stopCluster(cl)
  Reduce(rbind, xx)
}
