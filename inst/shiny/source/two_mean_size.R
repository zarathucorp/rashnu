#' Sample Size or Power for Two-Sample Mean Test
#'
#' Calculates sample size or power for a two-sample mean test.
#'
#' @param muA Numeric. True mean of group A.
#' @param muB Numeric. True mean of group B.
#' @param delta Numeric (optional). Margin for `"non-inferiority"` or `"equivalence test"`. Required for `"non-inferiority"` or `"equivalence"` test.
#' @param kappa Numeric. Ratio of sample sizes (nA/nB). Default is 1.
#' @param sd Numeric (optional). Standard deviation. Required for `"2-side"`, `"non-inferiority"` or `"equivalence"` test.
#' @param sdA Numeric (optional). Standard deviation of group A. Required for `"1-side"` test.
#' @param sdB Numeric (optional). Standard deviation of group B. Required for `"1-side"` test.
#' @param alpha Numeric. Type I error rate.
#' @param beta Numeric (optional). Type II error rate. Required for sample size calculation.
#' @param nA Integer (optional). Sample size for group A. Required for power calculation of `"1-side"` test.
#' @param nB Integer (optional). Sample size for group B. Required for power calculation of `"2-side"`, `"non-inferiority"` or `"equivalence"` test.
#' @param test_type Character. `"2-side"`, `"1-side"`, `"non-inferiority"`, or `"equivalence"`. Default is `"2-side"`.
#'
#' @return Numeric. Returns sample size (if `beta` is given), or power (if `nA`/`nB` is given).
#'
#' @note
#' Only one of `beta` (for sample size calculation) or `nA`/`nB` (for power calculation) should be specified.
#'
#' Required arguments by `test_type`:
#' - `"2-side"`:
#'   - For sample size: `muA`, `muB`, `sd`, `alpha`, `beta`
#'   - For power: `muA`, `muB`, `sd`, `alpha`, `nB`
#'
#' - `"1-side"`:
#'   - For sample size: `muA`, `muB`, `sdA`, `sdB`, `alpha`, `beta`
#'   - For power: `muA`, `muB`, `sdA`, `sdB`, `alpha`, `nA`
#'
#' - `"non-inferiority"`/`"equivalence"`:
#'   - For sample size: `muA`, `muB`, `delta`, `sd`, `alpha`, `beta`
#'   - For power: `muA`, `muB`, `delta`, `sd`, `alpha`, `nB`
#'
#' @examples
#' # Sample size for `"2-side"` test
#' two_mean_size(muA = 5, muB = 10, kappa = 1, sd = 10,
#'               alpha = 0.05, beta = 0.2, test_type = "2-side")
#'
#' # Power of `"2-side"` test
#' two_mean_size(muA = 5, muB = 10, kappa = 1, sd = 10,
#'               alpha = 0.05, nB = 63, test_type = "2-side")
#'
#' # Sample size for `"1-side"` test
#' two_mean_size(muA = 132.86, muB = 127.44, kappa = 2, sdA = 15.34, sdB = 18.23,
#'               alpha = 0.05, beta = 0.2, test_type = "1-side")
#'
#' # Power of `"1-sided"` test
#' two_mean_size(muA = 132.86, muB = 127.44, kappa = 2, sdA = 15.34, sdB = 18.23,
#'               alpha = 0.05, nA = 85, test_type = "1-side")
#'
#' # Sample size for `"non-inferiority"` test
#' two_mean_size(muA = 5, muB = 5, delta = 5, kappa = 1, sd = 10,
#'               alpha = 0.05, beta = 0.2, test_type = "non-inferiority")
#'
#' # Power of `"non-inferiority"` test
#' two_mean_size(muA = 5, muB = 5, delta = 5, kappa = 1, sd = 10,
#'               alpha = 0.05, nB = 50, test_type = "non-inferiority")
#'
#' # Sample size for `"equivalence"` test
#' two_mean_size(muA = 5, muB = 4, delta = 5, kappa = 1, sd = 10,
#'               alpha = 0.05, beta = 0.2, test_type = "equivalence")
#'
#' # Power of `"equivalence"` test
#' two_mean_size(muA = 5, muB = 4, delta = 5, kappa = 1, sd = 10,
#'               alpha = 0.05, nB = 108, test_type = "equivalence")
#'
#' @export
two_mean_size <- function(muA, muB, delta = NULL, kappa = 1, sd = NULL, sdA = NULL, sdB = NULL, alpha, beta = NULL, nA = NULL, nB = NULL, test_type = "2-side") {
  if (!is.null(beta)) {
    if (test_type == "2-side") {
      return(ceiling((1+1/kappa)*(sd*(qnorm(1-alpha/2)+qnorm(1-beta))/(muA-muB))^2))
    } else if (test_type == "1-side") {
      return(ceiling((sdA^2+sdB^2/kappa)*((qnorm(1-alpha)+qnorm(1-beta))/(muA-muB))^2))
    } else if (test_type == "non-inferiority") {
      return(ceiling((1+1/kappa)*(sd*(qnorm(1-alpha)+qnorm(1-beta))/(muA-muB-delta))^2))
    } else if (test_type == "equivalence") {
      return(ceiling((1+1/kappa)*(sd*(qnorm(1-alpha)+qnorm(1-beta/2))/(abs(muA-muB)-delta))^2))
    }
  } else if (!is.null(nA)) {
    if (test_type == "1-side") {
      return(pnorm((muA-muB)/sqrt(sdA^2+sdB^2/kappa)*sqrt(nA)-qnorm(1-alpha)))
    }
  } else if (!is.null(nB)) {
    if (test_type == "2-side") {
      return(pnorm((muA-muB)/(sd*sqrt((1+1/kappa)/nB))-qnorm(1-alpha/2))+pnorm(-(muA-muB)/(sd*sqrt((1+1/kappa)/nB))-qnorm(1-alpha/2)))
    } else if (test_type == "non-inferiority") {
      return(pnorm((muA-muB-delta)/(sd*sqrt((1+1/kappa)/nB))-qnorm(1-alpha))+pnorm(-(muA-muB-delta)/(sd*sqrt((1+1/kappa)/nB))-qnorm(1-alpha)))
    } else if (test_type == "equivalence") {
      return(2*(pnorm((abs(muA-muB)-delta)/(sd*sqrt((1+1/kappa)/nB))-qnorm(1-alpha))+pnorm(-(abs(muA-muB)-delta)/(sd*sqrt((1+1/kappa)/nB))-qnorm(1-alpha)))-1)
    }
  }
}
