% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Visuals.R
\name{plotClusters}
\alias{plotClusters}
\title{plot of clusters with weighted edges}
\usage{
plotClusters(
  clus,
  myModel,
  eWidth = 2,
  vSize = 30,
  vFontSize = 2,
  vColor = "cyan"
)
}
\arguments{
\item{clus}{clusters obtained from makeSNPClusters()}

\item{myModel}{cluster analyis model created using createModel()}

\item{eWidth}{width of edges}

\item{vSize}{size of vertices}

\item{vFontSize}{font size of vertex labels}

\item{vColor}{vertex colour}
}
\value{
graph object
}
\description{
plot of clusters with weighted edges
}
