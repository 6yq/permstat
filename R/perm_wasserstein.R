#' Wasserstein Distance Permutation Test
#'
#' Performs a permutation test based on the Wasserstein-1 distance between two numeric samples.
#' This test assesses whether two distributions differ in location or shape.
#'
#' @importFrom transport wasserstein1d
#' @export
#'
#' @param x Numeric vector, first sample.
#' @param y Numeric vector, second sample.
#' @param N Integer (default 10000), number of permutations.
#' @param seed Optional integer. If provided, sets random seed for reproducibility.
#'
#' @return A list with:
#'   \item{observed_stat}{Observed Wasserstein distance}
#'   \item{p_value}{Permutation p-value}
#'   \item{perm_distribution}{Vector of permuted distances}
#'
#' @examples
#' \dontrun{
#'   library(transport)
#'   perm_wasserstein(rnorm(20), rnorm(30, mean = 1), seed = 123)
#' }

perm_wasserstein <- function(x, y, N = 10000, seed = NULL) {
  stopifnot(is.numeric(x), is.numeric(y))
  if (!requireNamespace("transport", quietly = TRUE)) {
    stop("Package 'transport' is required. Please install it with install.packages('transport')")
  }
  if (!is.null(seed)) set.seed(seed)

  combined <- c(x, y)
  n <- length(x)

  stat_fun <- function(a, b) {
    transport::wasserstein1d(a, b)
  }

  T_obs <- stat_fun(x, y)

  perm_stats <- replicate(N, {
    idx <- sample(length(combined), n)
    stat_fun(combined[idx], combined[-idx])
  })

  p_val <- mean(perm_stats >= T_obs)

  list(
    observed_stat = T_obs,
    p_value = p_val,
    perm_distribution = perm_stats
  )
}
