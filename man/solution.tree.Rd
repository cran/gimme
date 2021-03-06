% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/solution.tree.R
\name{solution.tree}
\alias{solution.tree}
\title{Solution trees for multiple solutions gimme.}
\usage{
solution.tree(x,
              level     =  c("group", "individual"),
              cols      =  NULL,
              ids       =  "all",
              plot.tree =  FALSE)
}
\arguments{
\item{x}{A fitted gimme object.}

\item{level}{A character vector indicating what levels of the solution tree 
you would like returned.  Options are "group", "individual", or c("group", "individual").
Defaults to c("group", "individual").}

\item{cols}{A character vector indicating additional information to include in tree plot.
Options include "stage", "pruned", "rmsea", "nnfi", "cfi","srmr", "grp_sol",
"bic", "aic", "modularity." Defaults to NULL.}

\item{ids}{A character vector indicating the names of subjects to print.  Defaults to "all."}

\item{plot.tree}{Logical.  If TRUE, plot of tree is produced.  Defaults to FALSE.}
}
\description{
This function allows for the exploration of divergences in multiple
solutions gimme for both the group and individuals.
}
\details{
solution.tree
}
