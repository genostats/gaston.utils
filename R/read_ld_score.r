read.ld.score  <- function(filename, chr, pos) {
  if(missing(chr) | missing(pos)) 
    D <- .Call("read_ld_score_all", PACKAGE = "gaston.utils", filename)
  else
    D <- .Call("read_ld_score_filt", PACKAGE = "gaston.utils", filename, chr, pos)

  data.frame(D)
}
