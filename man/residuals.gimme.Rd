% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/residuals.gimme.R
\name{residuals.gimme}
\alias{residuals.gimme}
\title{GIMME Residuals.}
\usage{
residuals.gimme(x)
}
\arguments{
\item{x}{A fitted gimme object.}
}
\value{
List of two lists of data frames. \describe{ \item{residuals}{List of
  the unstandardized residuals per subject.}
  \item{standardized.residuals}{List of the standardized residuals per
  subject.}}
}
\description{
This function calculates the unstandardized and standardized
  residuals of a fitted gimme model.
}
\examples{
 \dontrun{
paths <- 'V2 ~ V1
          V3 ~ V4lag'

fit <- gimmeSEM(data     = simData,
                out      = "C:/simData_out",
                subgroup = TRUE,
                paths    = paths)

residuals <- residuals.gimme(fit)
 }
}
\author{
Sebastian Castro-Alvarez
}
\keyword{internal}
