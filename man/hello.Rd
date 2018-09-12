\name{geneFiltering}
\alias{geneFiltering}
\title{Gene Filtering for Pathway Analysis}
\usage{
simPathwaySystem(alpha, Rsq, nsample, rind, nx, nxPrime, ns, ntranscriptNoise, ng, ngPrime, nh, ngeneNoise, noiseSD)
}
\description{
Performs simulations for two data types and the corresponding continuous outcome based off of a multivariate normal distribution and latent variables.
}
\examples{
simResults <- simPathwaySystem(
  alpha= 0.35
  ,Rsq= 0.35
  ,nsample= 200
  ,rind= 0.35
  ,nx= 15
  ,nxPrime= 100
  ,ns= 50
  ,ntranscriptNoise= 100
  ,ng= 10
  ,ngPrime= 100
  ,nh= 30
  ,ngeneNoise= 100
  ,noiseSD= 0.5
)
}
