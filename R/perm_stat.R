#' Tail-Sensitive Permutation Test
#'
#' Performs a permutation test comparing the right-tail proportions of two numeric vectors.
#' The test statistic is the difference in proportions of values exceeding a specified quantile threshold.
#'
#' @importFrom stats quantile
#' @export
#'
#' @param x Numeric vector, first sample.
#' @param y Numeric vector, second sample.
#' @param q Numeric scalar between 0 and 1 (default 0.9), quantile threshold for defining "tail".
#' @param N Integer (default 10000), number of permutations.
#' @param seed Optional integer. If provided, sets random seed for reproducibility.
#'
#' @return A list with:
#'   \item{observed_stat}{Observed difference in right-tail proportions}
#'   \item{p_value}{Two-sided permutation p-value}
#'   \item{quantile_threshold}{The combined-sample quantile threshold used}
#'   \item{perm_distribution}{Vector of permuted statistics}
#'
#' @examples
#' x <- rlnorm(80)
#' y <- rlnorm(100)
#' perm_tail(x, y, seed = 42)

perm_tail <- function(x, y, q = 0.9, N = 10000, seed = NULL) {
  stopifnot(is.numeric(x), is.numeric(y))
  if (!is.null(seed)) set.seed(seed)

  combined <- c(x, y)
  n <- length(x)
  q_thresh <- quantile(combined, q)

  tail_stat <- function(a, b, qcut) {
    mean(a > qcut) - mean(b > qcut)
  }

  T_obs <- tail_stat(x, y, q_thresh)

  perm_stats <- replicate(N, {
    idx <- sample(length(combined), n)
    x_perm <- combined[idx]
    y_perm <- combined[-idx]
    tail_stat(x_perm, y_perm, q_thresh)
  })

  p_val <- mean(abs(perm_stats) >= abs(T_obs))

  list(
    observed_stat = T_obs,
    p_value = p_val,
    quantile_threshold = q_thresh,
    perm_distribution = perm_stats
  )
}


#' Median Difference Permutation Test
#'
#' Performs a two-sided permutation test comparing the medians of two numeric samples.
#' The test statistic is the difference in sample medians.
#'
#' @export
#'
#' @param x Numeric vector, first sample.
#' @param y Numeric vector, second sample.
#' @param N Integer (default 10000), number of permutations.
#' @param seed Optional integer. If provided, sets random seed for reproducibility.
#'
#' @return A list with:
#'   \item{observed_stat}{Observed difference in medians}
#'   \item{p_value}{Two-sided permutation p-value}
#'   \item{perm_distribution}{Vector of permuted statistics}
#'
#' @examples
#' perm_median(rnorm(20), rnorm(30), seed = 42)

perm_median <- function(x, y, N = 10000, seed = NULL) {
  stopifnot(is.numeric(x), is.numeric(y))
  if (!is.null(seed)) set.seed(seed)
  
  combined <- c(x, y)
  n <- length(x)

  median_stat <- function(a, b) {
    median(a) - median(b)
  }

  T_obs <- median_stat(x, y)

  perm_stats <- replicate(N, {
    idx <- sample(length(combined), n)
    x_perm <- combined[idx]
    y_perm <- combined[-idx]
    median_stat(x_perm, y_perm)
  })

  p_val <- mean(abs(perm_stats) >= abs(T_obs))

  list(
    observed_stat = T_obs,
    p_value = p_val,
    perm_distribution = perm_stats
  )
}
