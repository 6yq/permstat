% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/perm_wasserstein.R
\name{perm_wasserstein}
\alias{perm_wasserstein}
\title{Wasserstein Distance Permutation Test}
\usage{
perm_wasserstein(x, y, N = 10000, seed = NULL)
}
\arguments{
\item{x}{Numeric vector, first sample.}

\item{y}{Numeric vector, second sample.}

\item{N}{Integer (default 10000), number of permutations.}

\item{seed}{Optional integer. If provided, sets random seed for reproducibility.}
}
\value{
A list with:
\item{observed_stat}{Observed Wasserstein distance}
\item{p_value}{Permutation p-value}
\item{perm_distribution}{Vector of permuted distances}
}
\description{
Performs a permutation test based on the Wasserstein-1 distance between two numeric samples.
This test assesses whether two distributions differ in location or shape.
}
\examples{
\dontrun{
  library(transport)
  perm_wasserstein(rnorm(20), rnorm(30, mean = 1), seed = 123)
}
}
