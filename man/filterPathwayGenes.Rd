\name{geneFiltering}
\alias{filterPathwayGenes}
\title{Gene Filtering Algorithm}
\usage{
filterPathwayGenes(x, g, y, numClusters_1, numClusters_2)
}
\description{
Filtering algorithm for two datatypes and a corresponding outcome.
}
\examples{
filterResults <- filterPathwayGenes(
  x=x
  ,g=g
  ,y=y
  ,numClusters_1=3
  ,numClusters_2=3
)
}
