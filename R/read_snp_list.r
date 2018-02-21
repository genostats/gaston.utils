read.snp.list  <- function(filename, id) {
  if(missing(id))
    D <- .Call("read_snp_list_all", PACKAGE = "gaston.utils", filename)
  else 
    D <- .Call("read_snp_list_filt_id", PACKAGE = "gaston.utils", filename, id)

  data.frame(D)
}
