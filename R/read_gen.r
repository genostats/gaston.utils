# this is just a test for reading a .gen file
read.gen.file <- function(x) {
  x <- path.expand(x)
  .Call('zz_read_gen_file', PACKAGE = "gaston.utils", x)
}

# cette fonction vérifie le format de tout le fichier
dim.gen.file <- function(x) {
  x <- path.expand(x)
  .Call('zz_gen_file_dim', PACKAGE = "gaston.utils", x)
}

# celle-ci ne lit que la première ligne
nb.inds.gen.file <- function(x) {
  x <- path.expand(x)
  .Call('zz_nb_inds_gen_file', PACKAGE = "gaston.utils", x)
}
