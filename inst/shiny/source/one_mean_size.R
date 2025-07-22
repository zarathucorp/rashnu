#' Sample Size or Power for One-Sample Mean Test
#'
#' Calculates sample size or power for a one-sample mean test.
#'
#' @param mu Numeric. True mean.
#' @param mu0 Numeric. Null hypothesis mean.
#' @param delta Numeric (optional). Margin for `"non-inferiority"` or `"equivalence"` test. Required for `"non-inferiority"` or `"equivalence"` test.
#' @param sd Numeric. Standard deviation.
#' @param alpha Numeric. Type I error rate.
#' @param beta Numeric (optional). Type II error rate. Required for sample size calculation.
#' @param n Integer (optional). Sample size. Required for power calculation.
#' @param test_type Character. `"2-side"`, `"1-side"`, `"non-inferiority"`, or `"equivalence"`. Default is `"2-side"`.
#'
#' @return Numeric. Returns sample size (if `beta` is given), or power (if `n` is given).
#'
#' @note
#' Only one of `beta` (for sample size calculation) or `n` (for power calculation) should be specified.
#'
#' Required arguments by `test_type`:
#' - `"2-side"` / `"1-side"`:
#'   - For sample size: `mu`, `mu0`, `sd`, `alpha`, `beta`
#'   - For power: `mu`, `mu0`, `sd`, `alpha`, `n`
#'
#' - `"non-inferiority"` / `"equivalence"`:
#'   - For sample size: `mu`, `mu0`, `delta`, `sd`, `alpha`, `beta`
#'   - For power: `mu`, `mu0`, `delta`, `sd`, `alpha`, `n`
#'
#' @examples
#' # Sample size for `"2-side"` test
#' one_mean_size(mu = 2, mu0 = 1.5, sd = 1,
#'               alpha = 0.05, beta = 0.2, test_type = "2-side")
#'
#' # Power of `"2-side"` test
#' one_mean_size(mu = 2, mu0 = 1.5, sd = 1,
#'               alpha = 0.05, n = 32, test_type = "2-side")
#'
#' # Sample size for `"1-side"` test
#' one_mean_size(mu = 115, mu0 = 120, sd = 24,
#'               alpha = 0.05, beta = 0.2, test_type = "1-side")
#'
#' # Power of `"1-side"` test
#' one_mean_size(mu = 115, mu0 = 120, sd = 24,
#'               alpha = 0.05, n = 143, test_type = "1-side")
#'
#' # Sample size for `"non-inferiority"` test
#' one_mean_size(mu = 2, mu0 = 1.5, delta = -0.5, sd = 1,
#'               alpha = 0.05, beta = 0.2, test_type = "non-inferiority")
#'
#' # Power of `"non-inferiority"` test
#' one_mean_size(mu = 2, mu0 = 1.5, delta = -0.5, sd = 1,
#'               alpha = 0.05, n = 7, test_type = "non-inferiority")
#'
#' # Sample size for `"equivalence"` test
#' one_mean_size(mu = 2, mu0 = 2, delta = 0.05, sd = 0.1,
#'               alpha = 0.05, beta = 0.2, test_type = "equivalence")
#'
#' # Power of `"equivalence"` test
#' one_mean_size(mu = 2, mu0 = 2, delta = 0.05, sd = 0.1,
#'               alpha = 0.05, n = 35, test_type = "equivalence")
#'
#' @export
one_mean_size <- function(mu, mu0, delta = NULL, sd, alpha, beta = NULL, n = NULL, test_type = "2-side") {
  if (!is.null(beta)) {
    if (test_type == "2-side") {
      return(ceiling((sd*(qnorm(1-alpha/2)+qnorm(1-beta))/(mu-mu0))^2))
    } else if (test_type == "1-side") {
      return(ceiling((sd*(qnorm(1-alpha)+qnorm(1-beta))/(mu-mu0))^2))
    } else if (test_type == "non-inferiority") {
      return(ceiling((sd*(qnorm(1-alpha)+qnorm(1-beta))/(mu-mu0-delta))^2))
    } else if (test_type == "equivalence") {
      return(ceiling((sd*(qnorm(1-alpha)+qnorm(1-beta/2))/(delta-abs(mu-mu0)))^2))
    }
  } else if (!is.null(n)) {
    if (test_type == "2-side") {
      return(pnorm((mu-mu0)/sd*sqrt(n)-qnorm(1-alpha/2))+pnorm(-(mu-mu0)/sd*sqrt(n)-qnorm(1-alpha/2)))
    } else if (test_type == "1-side") {
      return(pnorm(abs((mu-mu0)/sd*sqrt(n))-qnorm(1-alpha)))
    } else if (test_type == "non-inferiority") {
      return(pnorm((mu-mu0-delta)/sd*sqrt(n)-qnorm(1-alpha))+pnorm(-(mu-mu0-delta)/sd*sqrt(n)-qnorm(1-alpha)))
    } else if (test_type == "equivalence") {
      return(2*(pnorm((abs(mu-mu0)-delta)/sd*sqrt(n)-qnorm(1-alpha))+pnorm(-(abs(mu-mu0)-delta)/sd*sqrt(n)-qnorm(1-alpha)))-1)
    }
  }
}
