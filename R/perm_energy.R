#' Energy Distance Permutation Test
#'
#' Performs a permutation test using the energy distance between two numeric samples.
#' This test measures whether the two distributions are statistically distinguishable.
#'
#' @importFrom stats dist
#' @export
#'
#' @param x Numeric vector, first sample.
#' @param y Numeric vector, second sample.
#' @param N Integer (default 10000), number of permutations.
#' @param seed Optional integer. If provided, sets random seed for reproducibility.
#'
#' @return A list with:
#'   \item{observed_stat}{Observed energy distance}
#'   \item{p_value}{Permutation p-value}
#'   \item{perm_distribution}{Vector of permuted distances}
#'
#' @examples
#' perm_energy(rnorm(20), rnorm(30, mean = 1), seed = 123)

perm_energy <- function(x, y, N = 10000, seed = NULL) {
  stopifnot(is.numeric(x), is.numeric(y))
  if (!is.null(seed)) set.seed(seed)

  n <- length(x)
  m <- length(y)
  combined <- c(x, y)
  sizes <- c(n, m)

  # Pairwise Euclidean distance
  energy_stat <- function(a, b) {
    ab <- dist(c(a, b))
    ab <- as.matrix(ab)
    n <- length(a)
    m <- length(b)

    d_ab <- mean(ab[1:n, (n+1):(n+m)])
    d_aa <- mean(ab[1:n, 1:n][upper.tri(matrix(1, n, n))])
    d_bb <- mean(ab[(n+1):(n+m), (n+1):(n+m)][upper.tri(matrix(1, m, m))])
    2 * d_ab - d_aa - d_bb
  }

  T_obs <- energy_stat(x, y)

  perm_stats <- replicate(N, {
    idx <- sample(length(combined), n)
    energy_stat(combined[idx], combined[-idx])
  })

  p_val <- mean(perm_stats >= T_obs)

  list(
    observed_stat = T_obs,
    p_value = p_val,
    perm_distribution = perm_stats
  )
}
