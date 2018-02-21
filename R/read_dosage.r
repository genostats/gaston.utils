# this is just a test for reading a .dose file (or dosages in a VCF)
read.dosage.file <- function(filename) {
  filename <- path.expand(filename)
  .Call('zz_read_dose_file', PACKAGE = "gaston.utils", filename)
}

# cette fonction vérifie le format de tout le fichier
dim.dosage.file <- function(filename) {
  filename <- path.expand(filename)
  .Call('zz_dose_file_dim', PACKAGE = "gaston.utils", filename)
}

# celle-ci ne lit que la première ligne
nb.inds.dosage.file <- function(filename) {
  filename <- path.expand(filename)
  .Call('zz_nb_inds_dose_file', PACKAGE = "gaston.utils", filename)
}
