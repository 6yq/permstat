% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/perm_stat.R
\name{perm_tail}
\alias{perm_tail}
\title{Tail-Sensitive Permutation Test}
\usage{
perm_tail(x, y, q = 0.9, N = 10000, seed = NULL)
}
\arguments{
\item{x}{Numeric vector, first sample.}

\item{y}{Numeric vector, second sample.}

\item{q}{Numeric scalar between 0 and 1 (default 0.9), quantile threshold for defining "tail".}

\item{N}{Integer (default 10000), number of permutations.}

\item{seed}{Optional integer. If provided, sets random seed for reproducibility.}
}
\value{
A list with:
\item{observed_stat}{Observed difference in right-tail proportions}
\item{p_value}{Two-sided permutation p-value}
\item{quantile_threshold}{The combined-sample quantile threshold used}
\item{perm_distribution}{Vector of permuted statistics}
}
\description{
Performs a permutation test comparing the right-tail proportions of two numeric vectors.
The test statistic is the difference in proportions of values exceeding a specified quantile threshold.
}
\examples{
x <- rlnorm(80)
y <- rlnorm(100)
perm_tail(x, y, seed = 42)
}
