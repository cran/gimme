\name{gimme-package}
\alias{gimme-package}
\alias{gimme-package}
\docType{package}
\title{
gimme: group iterative multiple model estimation
}
\description{
This package contains functions to identify group- and individual-level structural equation models. 
}
\details{
\tabular{ll}{
Package: \tab gimme\cr
Type: \tab Package\cr
Version: \tab 0.1-1\cr
Date: \tab 2014-10-15\cr
License: \tab GPL-2 \cr
}

Researchers across varied domains gather multivariate data for each individual unit of study across multiple occasions of measurement. Generally referred to as time series (or in the social sciences, intensive longitudinal) data, examples include psychophysiological processes such as neuroimaging and heart rate variability, daily diary studies, and observational coding of social interactions among dyads. A primary goal for acquiring these data is to understand temporal processes. 

The gimme package contains several functions for use with these data. These functions include \code{\link{gimme}}, which provides both group- and individual-level results by looking across individuals for patterns of relations among variables. A function that provides group-level results, \code{\link{aggregate}}, is included, as well as a function that provides individual-level results, \code{\link{indSEM}}. Finally, \code{\link{semigimme}} allows the user to specify a priori paths with which to begin estimation. 

The major functions within the gimme package all require the user to specify the directory containing the data and a directory for output to be stored. 
}
\author{

Stephanie Lane  <slane@unc.edu>
}
\references{
Gates, K.M. & Molenaar, P.C.M. (2012). Group search algorithm recovers effective connectivity maps for individuals in homogeneous and heterogeneous samples. NeuroImage, 63, 310-319.
}

\keyword{ package }

