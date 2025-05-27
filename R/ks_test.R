#' Two-Sample Kolmogorov-Smirnov Test with Exact or Permutation p-Value
#'
#' Performs a two-sample Kolmogorov-Smirnov test comparing the empirical distributions
#' of two numeric vectors. If exact p-values are not available due to sample size, ties,
#' or one-sided alternatives, a permutation test is performed using the same KS statistic.
#'
#' @importFrom stats ks.test
#' @export
#'
#' @param x Numeric vector, first sample.
#' @param y Numeric vector, second sample.
#' @param alternative Character string specifying the alternative hypothesis:
#'   \code{"two.sided"} (default), \code{"less"}, or \code{"greater"}.
#' @param N Integer (default 10000), number of permutations if exact p-value is not available.
#' @param seed Optional integer seed for reproducibility of permutation sampling.
#'
#' @return A list with:
#'   \item{method}{Test method used: "KS test (exact)" or "KS test (permutation)"}
#'   \item{observed_stat}{Observed KS test statistic}
#'   \item{p_value}{Two-sided or one-sided p-value depending on the alternative}
#'   \item{type}{Type of p-value: "exact" or "permutation"}
#'   \item{perm_distribution}{(If permutation used) Vector of permuted test statistics}
#'
#' @examples
#' x <- rnorm(50)
#' y <- rnorm(40, mean = 0.3)
#' perm_ks(x, y, alternative = "greater")

perm_ks <- function(x, y, alternative = "two.sided", N = 10000, seed = NULL) {
  stopifnot(is.numeric(x), is.numeric(y))
  n <- length(x)
  m <- length(y)
  combined <- c(x, y)

  has_ties <- any(duplicated(combined))
  one_sided <- alternative %in% c("less", "greater")
  too_large <- n * m >= 10000

  use_exact <- !(has_ties || one_sided || too_large)

  # Function to extract KS statistic safely
  ks_stat <- function(a, b) {
    suppressWarnings(as.numeric(ks.test(a, b, alternative = alternative)$statistic))
  }

  if (use_exact) {
    # Exact test via ks.test
    ks <- suppressWarnings(ks.test(x, y, alternative = alternative, exact = TRUE))
    return(list(
      method = "KS test (exact)",
      observed_stat = as.numeric(ks$statistic),
      p_value = ks$p.value,
      type = "exact"
    ))
  } else {
    # Permutation test
    if (!is.null(seed)) set.seed(seed)
    T_obs <- ks_stat(x, y)
    perm_stats <- replicate(N, {
      idx <- sample(n + m)
      x_perm <- combined[idx[1:n]]
      y_perm <- combined[idx[(n + 1):(n + m)]]
      ks_stat(x_perm, y_perm)
    })

    p_val <- mean(perm_stats >= T_obs)

    return(list(
      method = "KS test (permutation)",
      observed_stat = T_obs,
      p_value = p_val,
      type = "permutation",
      perm_distribution = perm_stats
    ))
  }
}
