#!/usr/bin/env Rscript

# TODO: could allow install location to be supplied in call of script

##### Load Packages Needed #####
# includes dependencies as defined in utilities.R @import calls
my_packages <- c("glue", "devtools", "roxygen2", "usethis", "ggplot2", "magick",
                 "Seurat", "reticulate", "tidyverse")

# define repo to install from
options(repos = c(CRAN = "http://cran.rstudio.com"))

# install / import using CRAN or Bioconductor, as appropriate
for (pkg in my_packages) {
  if (!require(pkg, character.only = TRUE)) {
    tryCatch(
      {install.packages(pkg, lib = .libPaths()[1])
       invisible(lapply(pkg, library, character.only = TRUE))},
      error = function(e) {
        message("Installation via install.packages() failed. Trying
                 BiocManager.")
        BiocManager::install(pkg, lib = .libPaths()[1])
        invisible(lapply(pkg, library, character.only = TRUE))
      },
      warning = function(w) {
        message(w)
      }
    )
  }
}

##### Alter the following, if appropriate: #####

# RELATIVE PATH for github (call script anywhere within repo file structure):
# get location of base_dir irrespective of where script is called
initial_options <- commandArgs(trailingOnly = FALSE)
file_arg_name <- "--file="
script_name <- sub(file_arg_name, "",
                   initial_options[grep(file_arg_name, initial_options)])
# Get the directory of the script
parent_path <- dirname(normalizePath(script_name))

package_name <- "integration"

# annotated function file to include in package (include roxygen comments):
fxns_file <- "utilities.R"

##### Main #####
# overwrite the environment's challenge_nested_project function
challenge_nested_project <- function(path, name) {
  return() # return nothing (overwrite)
}

# Unlock the environment and bindings of the usethis package
rlang::env_unlock(env = asNamespace("usethis"))
rlang::env_binding_unlock(env = asNamespace("usethis"))

# Assign the challenge_nested_project function to the usethis package namespace
assign("challenge_nested_project", challenge_nested_project,
       envir = asNamespace("usethis"))

# Lock the bindings and environment of the usethis package again
rlang::env_binding_lock(env = asNamespace("usethis"))
rlang::env_lock(asNamespace("usethis"))

# Create project in <parent_path> folder
create_project(path = package_name, open = TRUE, rstudio = TRUE)

# ##### Alter the following, if appropriate: #####
# # Create DESCRIPTION file in <parent_path>/<package_name> folder
use_description(roxygen = TRUE,
                fields = list(Title =
                "Integrating, Merging, and Querying Datasets",
                Description = glue("Functions relevant for the integration,",
                " of datasets and using said integrated datasets as a",
                " reference from which cells can be mapped and queried for",
                " appropriate cell type labels. Includes plotting defaults,",
                " a means of overlaying ggplots, etc."),
                `Authors@R` = NULL,
                License = "Talk to Kimberly Aldinger @ SCRI",
                Authors = "Lindsey Sydnor @ lindsey.sydnor@gmail.com",
                Version = "1.0", Encoding = "UTF-8"))

# Copy contents of fxns_file (functions and roxygen annotations, to package dir)
# See "Template" at bottom of page for example of fxns_file contents
file.copy(file.path(parent_path, fxns_file),
          file.path(parent_path, package_name, "R", fxns_file),
          overwrite = TRUE)

# NOTE, will overwrite in package dir if you make edits.
# Install package to .libPaths()[1]
roxygenise()
install(file.path(parent_path, package_name))

# check to see if it's there
print("Contents of libPaths (check for integration installation):")
print(list.files(.libPaths()[1]))

##### Template #####

# #' save_fxn
# #'
# #' @description Save a ggplot object to the appropriate location(s)
# #' @param obj (ggplot obj): ggplot object to save
# #' @param file (string): filename to save ggplot as (include file extension:
# #'  'png' only)
# #' @param w (int): width of PNG to save in inches
# #' @param h (int): height of PNG to save in inches
# #' @param dirs (string or character vector): directory (or directories) to save
# #'  the plot in
# #' @import ggplot2 glue
# #' @return N/A
# #' @examples 
# #' save_fxn(main_plt, file = "tmp_plt1.png", dirs = c("~/endothelial/images",
# #' glue("/active/aldinger_k/kimslab_temp/scRNA-seq/brain-vasc/vascular-dev/"
# #' "lsyd/images")))
# #' @export
# save_fxn <- function(obj, file, w = 8, h = 5,
#                      dirs = c(glue("/active/aldinger_k/kimslab_temp/scRNA-seq/",
#                                    "brain-vasc/vascular-dev/lsyd/images/"),
#                               "~/endothelial/images/")) {
#   for (d in dirs) {
#     dir.create(d, showWarnings = FALSE, recursive = TRUE)
#     ggsave(filename = file.path(d, file), plot = obj, device = "png",
#            units = "in", width = w, height = h)
#   }
#   # save 1, copy over rest (faster)
#   if (length(dirs) > 1) {
#     for (i in 2:length(dirs)) {
#       file.copy(from = paste0(dirs[1], file),
#                 to = paste0(dirs[i], file), overwrite = TRUE)
#     }
#   }
# }
