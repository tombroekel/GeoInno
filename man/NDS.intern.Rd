% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/complexity_calculation.R
\name{NDS.intern}
\alias{NDS.intern}
\title{individual network diversity score (iNDS)}
\usage{
NDS.intern(node_x = x, g, reps_i)
}
\arguments{
\item{node_x}{The randomly sampled node (position index in igraph object) for which the partial network is to be extracted by means of a random walk.}

\item{g}{The igraph object. Usually, a binarized version of the combinatorial network of CPC classes co-occurring on patents.}

\item{reps_i}{The number of repetitions used in the bootstrap, default set to 200.}
}
\value{
Returns the value of the iNDS measure.
}
\description{
NDS.intern() calculates the individual network diversity score as defined by \insertCite{Emmert-Streib2012;textual}{GeoInno}. It is in the calculation of the structural diversity complexity measure of \insertCite{Broekel2019;textual}{GeoInno}.
}
\examples{
my.graph <- igraph::random.graph.game(p.or.m = 1/10, n=10)
NDS.intern(node_x=1, g = my.graph, reps_i=10)
}
