
fit_for_fst <- function(x, X) {
  if(!is(x, "bed.matrix")) 
    stop("...")

  .Call("gg_fit_allelic_freq", PACKAGE = "gaston.utils", x@bed, x@p, X, 0, ncol(x)-1, 1e-4);
}

Fst <- function(x, X) {
  if(!is(x, "bed.matrix"))
    stop("...")

  .Call("gg_fst", PACKAGE = "gaston.utils", x@bed, x@p, X, 0, ncol(x)-1, 1e-4);
}

