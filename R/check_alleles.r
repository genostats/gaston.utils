check.alleles <- function(A1, A2, B1, B2) {
  data.frame( .Call("gg_check_alleles", PACKAGE="gaston.utils", A1, A2, B1, B2) )
}

