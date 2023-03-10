% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/complexity_calculation.R
\name{NDS.intern}
\alias{NDS.intern}
\title{individual network diversity score (iNDS)}
\usage{
NDS.intern(x, g, reps, ...)
}
\arguments{
\item{g}{The igraph object. Usually, a binarized version of the combinatorial network of CPC classes co-occurring on patents.}

\item{reps}{The number of repetitions used in the bootstrap, default set to 200.}

\item{s}{The (random) sample of nodes (position indizes in igraph object) for which the partial networks are to be extracted by means of a random walk.}

\item{node.sample}{The number of nodes sampled in the Network Diversity Score calculation, set to 125, see \insertCite{Emmert-Streib2012;textual}{GeoInno}}
}
\value{
Returns the value of the iNDS measure.
}
\description{
NDS.intern() calculates the individual network diversity score as defined by \insertCite{Emmert-Streib2012;textual}{GeoInno}. It is in the calculation of the structural diversity complexity measure of \insertCite{Broekel2019;textual}{GeoInno}.
}
\examples{
my.graph <- igraph::random.graph.game(p.or.m = 1/10, n=10)
NDS.intern(g = my.graph, s=c(1:10), node.sample=10, reps=10)
}
