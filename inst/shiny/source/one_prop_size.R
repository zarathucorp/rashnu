#' Sample Size or Power for One-Sample Proportion Test
#'
#' Calculates sample size or power for a one-sample proportion test.
#'
#' @param p Numeric. True proportion.
#' @param p0 Numeric. Null hypothesis proportion.
#' @param delta Numeric (optional). Margin for `"non-inferiority"` or `"equivalence"` test. Required for `"non-inferiority"` or `"equivalence"` test.
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
#' - `"2-side"`/`"1-side`:
#'   - For sample size: `p`, `p0`, `alpha`, `beta`
#'   - For power: `p`, `p0`, `alpha`, `n`
#'
#' - `"non-inferiority"`/`"equivalence"`:
#'   - For sample size: `p`, `p0`, `delta`, `alpha`, `beta`
#'   - For power: `p`, `p0`, `sdA`, `delta`, `alpha`, `n`
#'
#' @examples
#' # Sample size for `"2-side"` test
#' one_prop_size(p = 0.5, p0 = 0.3,
#'               alpha = 0.05, beta = 0.2, test_type = "2-side")
#'
#' # Power of `"2-side"` test
#' one_prop_size(p = 0.5, p0 = 0.3,
#'               alpha = 0.05, n = 50, test_type = "2-side")
#'
#' # Sample size for `"1-side"` test
#' one_prop_size(p = 0.05, p0 = 0.02,
#'               alpha = 0.05, beta = 0.2, test_type = "1-side")
#'
#' # Power of `"1-sided"` test
#' one_prop_size(p = 0.05, p0 = 0.02,
#'               alpha = 0.05, n = 191, test_type = "1-side")
#'
#' # Sample size for `"non-inferiority"` test
#' one_prop_size(p = 0.5, p0 = 0.3, delta = -0.1,
#'               alpha = 0.05, beta = 0.2, test_type = "non-inferiority")
#'
#' # Power of `"non-inferiority"` test
#' one_prop_size(p = 0.5, p0 = 0.3, delta = -0.1,
#'               alpha = 0.05, n = 18, test_type = "non-inferiority")
#'
#' # Sample size for `"equivalence"` test
#' one_prop_size(p = 0.6, p0 = 0.6, delta = 0.2,
#'               alpha = 0.05, beta = 0.2, test_type = "equivalence")
#'
#' # Power of `"equivalence"` test
#' one_prop_size(p = 0.6, p0 = 0.6, delta = 0.2,
#'               alpha = 0.05, n = 52, test_type = "equivalence")
#'
#' @export
one_prop_size <- function(p, p0, delta = NULL, alpha, beta = NULL, n = NULL, test_type = "2-side") {
  if (!is.null(beta)) {
    if (test_type == "2-side") {
      return(ceiling(p*(1-p)*((qnorm(1-alpha/2)+qnorm(1-beta))/(p-p0))^2))
    } else if (test_type == "1-side") {
      return(ceiling(p0*(1-p0)*((qnorm(1-alpha)+qnorm(1-beta)*sqrt(p*(1-p)/p0/(1-p0)))/(p-p0))^2))
    } else if (test_type == "non-inferiority") {
      return(ceiling(p*(1-p)*((qnorm(1-alpha)+qnorm(1-beta))/(p-p0-delta))^2))
    } else if (test_type == "equivalence") {
      return(ceiling(p*(1-p)*((qnorm(1-alpha)+qnorm(1-beta/2))/(abs(p-p0)-delta))^2))
    }
  } else if (!is.null(n)) {
    if (test_type == "2-side") {
      return(pnorm((p-p0)/sqrt(p*(1-p)/n)-qnorm(1-alpha/2))+pnorm(-(p-p0)/sqrt(p*(1-p)/n)-qnorm(1-alpha/2)))
    } else if (test_type == "1-side") {
      return(pnorm(sqrt(p0*(1-p0)/p/(1-p))*(abs((p-p0)/sqrt(p0*(1-p0)/n))-qnorm(1-alpha))))
    } else if (test_type == "non-inferiority") {
      return(pnorm((p-p0-delta)/sqrt(p*(1-p)/n)-qnorm(1-alpha))+pnorm(-(p-p0-delta)/sqrt(p*(1-p)/n)-qnorm(1-alpha)))
    } else if (test_type == "equivalence") {
      return(2*(pnorm((abs(p-p0)-delta)/sqrt(p*(1-p)/n)-qnorm(1-alpha))+pnorm(-(abs(p-p0)-delta)/sqrt(p*(1-p)/n)-qnorm(1-alpha)))-1)
    }
  }
}
