% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MIXTURE.DEBUG_V0.1.R
\name{GetMixture}
\alias{GetMixture}
\title{GetMixture
returns the regresssion coefficient matrix B as absolute values or as proportion (normalized) values}
\usage{
GetMixture(obj, type = c("proportion", "absolute"))
}
\arguments{
\item{obj}{an object of class MIXTURE (see MIXTURE)}

\item{type}{a character string "proportion" or "absolute". short word is allowed}
}
\value{
a SxK data matrix with K cell types (according to the signature matrix) and S subjects
}
\description{
GetMixture
returns the regresssion coefficient matrix B as absolute values or as proportion (normalized) values
}
