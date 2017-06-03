# this is just a test of the gzstream library
gzcat <- function(x) {
  x <- path.expand(x)
  .Call('zz_gzcat', PACKAGE = "gaston.utils", x)
}


