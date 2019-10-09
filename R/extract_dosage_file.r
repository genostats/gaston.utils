
extract.dosage.file <- function(filename, outname, keep) {
  filename <- path.expand(filename)

  if( substr(outname, nchar(outname)-2, nchar(outname)) != ".gz" )
    outname <- paste0(outname, ".gz")
  outname <- path.expand(outname)

  .Call('extract_dosage_file', PACKAGE = "gaston.utils", filename, outname, keep)
}
