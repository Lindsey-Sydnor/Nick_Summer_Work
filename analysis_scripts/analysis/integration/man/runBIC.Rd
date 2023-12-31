% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{runBIC}
\alias{runBIC}
\title{runBIC}
\usage{
runBIC(
  obj,
  ...,
  nfeature_min = 1500,
  nfeature_max = 5000,
  increment = 500,
  regress = NULL,
  plt = TRUE,
  npcs = 50,
  file_path = "~/BIC.png"
)
}
\arguments{
\item{obj}{(Seurat object): \cr
Seurat object on which to run BIC analysis.}

\item{nfeature_min}{(int or double = 1500): \cr
Minimum number of variable features to test.}

\item{nfeature_max}{(int or double = 5000): \cr
Maximum number of variable features to test.}

\item{increment}{(int or double = 500): \cr
Test ideal VariableFeature nfeature number between nfeature_min and
nfeature_max by increments of \code{increment}}

\item{regress}{(NULL, string, or char vector): \cr
Vector of variables to regress in \link[Seurat]{ScaleData} call, if any.}

\item{file_path}{(string = "~/BIC.png"): \cr
Path to save BIC comparison png}
}
\value{
(named list) (optimal_nfeatures, bic) \cr
optimal_nfeatures: Ideal number of variable features as calculated via BIC.
\cr bic: vector of bic values for each number of variable features tested.
}
\description{
Find the ideal number of variable features to use for a Seurat object's
downstream analysis. \cr
Process: Run \link[Seurat]{FindVariableFeatures} to find \code{nfeature_max}
number of variable features. Then for each number of variable features:
\link[Seurat]{ScaleData}, run \link[Seurat]{SCTransform} on the variable
features, run PCA, find the proportion of explained variance of the total
dataset's variance for the first \code{npcs} components, calculate total variance
of the dataset's "scale.data" slot, and then calculate the BIC value for
that total PCA model. Compare BIC values across each number of variable
features. The number of variable features corresponding to the lowest BIC
value resulted in the "best" BIC model. If a plateau of BIC values is
reached, the lowest number of variable features will be preferred. A plot
of bic value per number of variable features will be produced in
\code{file_path}.
}
\examples{
 nvar_features <- runBIC(seurat_obj, file_path = glue("/active/aldinger_k/",
                         "kimslab_temp/scRNA-seq/brain-vasc/vascular-dev/",
                         "lsyd/test_save.png")
}
