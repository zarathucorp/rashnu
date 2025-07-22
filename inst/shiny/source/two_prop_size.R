#' Sample Size or Power for Two-Sample Proportion Test
#'
#' Calculates sample size or power for a two-sample proportion test.
#'
#' @param pA Numeric. True proportion of group A.
#' @param pB Numeric. True proportion of group B.
#' @param delta Numeric (optional). Margin for `"non-inferiority"` or `"equivalence"` test. Required for `"non-inferiority"` or `"equivalence"` test.
#' @param kappa Numeric. Ratio of sample sizes (nA/nB). Default is 1.
#' @param alpha Numeric. Type I error rate.
#' @param beta Numeric (optional). Type II error rate. Required for sample size calculation.
#' @param nB Integer (optional). Sample size for group B. Required for power calculation.
#' @param test_type Character. `"2-side"`, `"1-side"`, `"non-inferiority"`, or `"equivalence"`. Default is `"2-side"`.
#'
#' @return Numeric. Returns sample size (if `beta` is given), or power (if `nB` is given).
#'
#' @note
#' Only one of `beta` (for sample size calculation) or `nA`/`nB` (for power calculation) should be specified.
#'
#' Required arguments by `test_type`:
#' - `"2-side"`/`"1-side"`:
#'   - For sample size: `pA`, `pB`, `alpha`, `beta`
#'   - For power: `pA`, `pB`, `alpha`, `nB`
#'
#' - `"non-inferiority"`/`"equivalence"`:
#'   - For sample size: `pA`, `pB`, `delta`,  `alpha`, `beta`
#'   - For power: `pA`, `pB`, `delta`, `alpha`, `nB`
#'
#' @examples
#' # Sample size for `"2-side"` test
#' two_prop_size(pA = 0.65, pB = 0.85, kappa = 1,
#'               alpha = 0.05, beta = 0.2, test_type = "2-side")
#'
#' # Power of `"2-side"` test
#' two_prop_size(pA = 0.65, pB = 0.85, kappa = 1,
#'               alpha = 0.05, nB = 70, test_type = "2-side")
#'
#' # Sample size for `"1-side"` test
#' two_prop_size(pA = 0.65, pB = 0.85, kappa = 1,
#'               alpha = 0.05, beta = 0.2, test_type = "1-side")
#'
#' # Power of `"1-sided"` test
#' two_prop_size(pA = 0.65, pB = 0.85, kappa = 1,
#'               alpha = 0.05, nB = 55, test_type = "1-side")
#'
#' # Sample size for `"non-inferiority"` test
#' two_prop_size(pA = 0.85, pB = 0.65, delta = -0.1, kappa = 1,
#'               alpha = 0.05, beta = 0.2, test_type = "non-inferiority")
#'
#' # Power of `"non-inferiority"` test
#' two_prop_size(pA = 0.85, pB = 0.65, delta = -0.1, kappa = 1,
#'               alpha = 0.05, nB = 25, test_type = "non-inferiority")
#'
#' # Sample size for `"equivalence"` test
#' two_prop_size(pA = 0.65, pB = 0.85, delta = 0.05, kappa = 1,
#'               alpha = 0.05, beta = 0.2, test_type = "equivalence")
#'
#' # Power of `"equivalence"` test
#' two_prop_size(pA = 0.65, pB = 0.85, delta = 0.05, kappa = 1,
#'               alpha = 0.05, nB = 136, test_type = "equivalence")
#'
#' @export
two_prop_size <- function(pA, pB, delta = NULL, kappa = 1, alpha, beta = NULL, nB = NULL, test_type = "2-side") {
  if (!is.null(beta)) {
    if (test_type == "2-side") {
      return(ceiling((pA*(1-pA)/kappa+pB*(1-pB))*((qnorm(1-alpha/2)+qnorm(1-beta))/(pA-pB))^2))
    } else if (test_type == "1-side") {
      return(ceiling((pA*(1-pA)/kappa+pB*(1-pB))*((qnorm(1-alpha)+qnorm(1-beta))/(pA-pB))^2))
    } else if (test_type == "non-inferiority") {
      return(ceiling((pA*(1-pA)/kappa+pB*(1-pB))*((qnorm(1-alpha)+qnorm(1-beta))/(pA-pB-delta))^2))
    } else if (test_type == "equivalence") {
      return(ceiling((pA*(1-pA)/kappa+pB*(1-pB))*((qnorm(1-alpha)+qnorm(1-beta/2))/(abs(pA-pB)-delta))^2))
    }
  } else if (!is.null(nB)) {
    if (test_type == "2-side") {
      return(pnorm((pA-pB)/sqrt(pA*(1-pA)/nB/kappa+pB*(1-pB)/nB)-qnorm(1-alpha/2))+pnorm(-(pA-pB)/sqrt(pA*(1-pA)/nB/kappa+pB*(1-pB)/nB)-qnorm(1-alpha/2)))
    } else if (test_type == "1-side") {
      return(pnorm(abs((pA-pB)/sqrt(pA*(1-pA)/nB/kappa+pB*(1-pB)/nB))-qnorm(1-alpha)))
    } else if (test_type == "non-inferiority") {
      return(pnorm((pA-pB-delta)/sqrt(pA*(1-pA)/nB/kappa+pB*(1-pB)/nB)-qnorm(1-alpha))+pnorm(-(pA-pB-delta)/sqrt(pA*(1-pA)/nB/kappa+pB*(1-pB)/nB)-qnorm(1-alpha)))
    } else if (test_type == "equivalence") {
      return(2*(pnorm((abs(pA-pB)-delta)/sqrt(pA*(1-pA)/nB/kappa+pB*(1-pB)/nB)-qnorm(1-alpha))+pnorm(-(abs(pA-pB)-delta)/sqrt(pA*(1-pA)/nB/kappa+pB*(1-pB)/nB)-qnorm(1-alpha)))-1)
    }
  }
}
