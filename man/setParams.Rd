% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Model.R
\name{setParams}
\alias{setParams}
\title{set parameters for existing cluster analysis model}
\usage{
setParams(thisModel, lambda = -1, beta = -1, percCutOff = -1)
}
\arguments{
\item{thisModel}{cluster analyis model created using createModel()}

\item{lambda}{(optional) clock rate}

\item{beta}{(optional) transmission rate}

\item{percCutOff}{(optional) cut-off percentage}
}
\value{
model object
}
\description{
set parameters for existing cluster analysis model
}
\examples{
myModel <- setParams(myModel, lambda=0.5)
}
