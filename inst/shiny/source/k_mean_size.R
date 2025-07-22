#' Sample Size or Power Calculation for K Means
#'
#' Calculates sample size or power for a multiple-sample mean test.
#'
#' @param muA Numeric. True mean of group A.
#' @param muB Numeric. True mean of group B.
#' @param kappa Numeric. Ratio of sample sizes (nA/nB). Default is 1.
#' @param sd Numeric (optional). Standard deviation. Required for `"2-side"` test.
#' @param sdA Numeric (optional). Standard deviation of group A. Required for `"1-side"` test.
#' @param sdB Numeric (optional). Standard deviation of group B. Required for `"1-side"` test.
#' @param tau Integer. Number of comparisons.
#' @param alpha Numeric. Type I error rate.
#' @param beta Numeric (optional). Type II error rate (1 - power). Required for sample size calculation.
#' @param n Integer (optional). Sample size. Required for power calculation of `"2-side"` test.
#' @param nA Integer (optional). Sample size of group A. Required for power calculation of `"1-side"` test.
#' @param test_type Character. `"2-side"` or `"1-side"`. Default is `"2-side"`
#'
#' @return Numeric. Returns sample size (if `beta` is given), or power (if `n`/`nA` is given).
#'
#' @note
#' Only one of `beta` (for sample size calculation) or `n`/`nA` (for power calculation) should be specified.
#'
#' Required arguments by `test_type`:
#' - `"2-side"`:
#'   - For sample size: `muA`, `muB`, `sd`, `tau`, `alpha`, `beta`
#'   - For power: `muA`, `muB`, `sd`, `tau`, `alpha`, `n`
#'
#' -`"1-side"`:
#'   - For sample size: `muA`, `muB`, `sdA`, `sdB`, `tau`, `alpha`, `beta`
#'   - For power: `muA`, `muB`, `sdA`, `sdB`, `tau`, `alpha`, `nA`
#'
#' @examples
#' # Sample size for `"2-side"` test
#' k_mean_size(muA = 5, muB = 10, sd = 10, tau = 1,
#'             alpha = 0.05, beta = 0.2, test_type = "2-side")
#'
#' # Power of `"2-side"` test
#' k_mean_size(muA = 5, muB = 10, sd = 10, tau = 1,
#'             alpha = 0.05, n = 63, test_type = "2-side")
#'
#' # Sample size for `"1-side"` test
#' k_mean_size(muA = 132.86, muB = 127.44, kappa = 2, sdA = 15.34, sdB = 18.23, tau = 1,
#'             alpha = 0.05, beta = 0.2, test_type = "1-side")
#'
#' # Power of `"1-side"` test
#' k_mean_size(muA = 132.86, muB = 127.44, kappa = 2, sdA = 15.34, sdB = 18.23, tau = 1,
#'             alpha = 0.05, nA = 85, test_type = "1-side")
#'
#' @export
k_mean_size <- function(muA, muB, kappa = 1, sd = NULL, sdA = NULL, sdB = NULL, tau = 1, alpha, beta = NULL, n = NULL, nA = NULL, test_type = "2-side") {
  if (!is.null(beta)) {
    if (test_type == "2-side") {
      return(ceiling(2*(sd*(qnorm(1-alpha/(2/tau))+qnorm(1-beta))/(muA-muB))^2))
    } else if (test_type == "1-side") {
      return(ceiling((sdA^2+sdB^2/kappa)*((qnorm(1-alpha/tau)+qnorm(1-beta))/(muA-muB))^2))
    }
  } else if (!is.null(n)) {
    if (test_type == "2-side") {
      return(pnorm((muA-muB)/(sd*sqrt(2/n))-qnorm(1-alpha/(2/tau)))+pnorm(-(muA-muB)/(sd*sqrt(2/n))-qnorm(1-alpha/(2/tau))))
    }
  } else if (!is.null(nA)) {
    if (test_type == "1-side") {
      return(pnorm((muA-muB)/sqrt(sdA^2+sdB^2/kappa)*sqrt(nA)-qnorm(1-alpha/tau)))
    }
  }
}
