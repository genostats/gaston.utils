read.vcf.header <- function(filename) {
  .Call("gg_read_vcf_head", PACKAGE = "gaston.utils", filename);
}

