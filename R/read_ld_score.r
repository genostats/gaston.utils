read.ld.score  <- function(filename, id, chr, pos) {
  if(missing(id) & (missing(chr) | missing(pos)))
    D <- .Call("read_ld_score_all", PACKAGE = "gaston.utils", filename)
  else if(missing(id) & !missing(chr) & !missing(pos))
    D <- .Call("read_ld_score_filt_chr_pos", PACKAGE = "gaston.utils", filename, chr, pos)
  else if(!missing(id))
    D <- .Call("read_ld_score_filt_id", PACKAGE = "gaston.utils", filename, id)

  total.nb.snps <- D$total_nb_snps
  D$total_nb_snps <- NULL
  list( total.nb.snps = total.nb.snps, ld.scores = data.frame(D, stringsAsFactors = FALSE) )
}
