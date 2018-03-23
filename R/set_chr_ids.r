set.chr.ids <- function(L) {
  if(missing(L))
    L <- list( x = getOption("gaston.chr.x",  23L), 
               y = getOption("gaston.chr.y",  24L), 
              mt = getOption("gaston.chr.mt", 26L) )
  .Call('gg_set_chr_ids', PACKAGE = "gaston.utils", L)
}
