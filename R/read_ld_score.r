read.ld.score  <- function(filename, chr, pos) {

  D <- .Call("read_all", PACKAGE = "gaston.utils", filename)

  D
}
