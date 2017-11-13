# this is just a test for reading a .dose file (or dosages in a VCF)
read.dosage.file <- function(x) {
  x <- path.expand(x)
  .Call('zz_read_dose_file', PACKAGE = "gaston.utils", x)
}

# cette fonction vérifie le format de tout le fichier
dim.dosage.file <- function(x) {
  x <- path.expand(x)
  .Call('zz_dose_file_dim', PACKAGE = "gaston.utils", x)
}

# celle-ci ne lit que la première ligne
nb.inds.dosage.file <- function(x) {
  x <- path.expand(x)
  .Call('zz_nb_inds_dose_file', PACKAGE = "gaston.utils", x)
}
