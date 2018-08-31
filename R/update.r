#' Update the patchseqtools library
#' 
update_patchseqtools <- function() {
  devtools::install_github("AllenInstitute/patchseqtools", 
    auth_token = "802976690281f1483c40de46d0a07e9d01a3de08")
}
