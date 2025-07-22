#' Sample Size or Power for Odds Ratio Test
#'
#' Calculates sample size or power for odds ratio test.
#'
#' @param pA Numeric. True proportion of group A.
#' @param pB Numeric. True proportion of group B.
#' @param delta Numeric (optional). Margin for `"non-inferiority"` or `"equivalence"` test. Required for `"non-inferiority"` or `"equivalence"` test.
#' @param kappa Numeric. Ratio of sample sizes (nA/nB). Default is 1.
#' @param alpha Numeric. Type I error rate.
#' @param beta Numeric (optional). Type II error rate. Required for sample size calculation.
#' @param nB Integer (optional). Sample size for group B. Required for power calculation.
#' @param test_type Character. `"equality"`, `"non-inferiority"`, or `"equivalence"`. Default is `"2-side"`.
#'
#' @return Numeric. Returns sample size (if `beta` is given), or power (if `nB` is given).
#'
#' @note
#' Only one of `beta` (for sample size calculation) or `nB` (for power calculation) should be specified.
#'
#' Required arguments by `test_type`:
#' - `"equality"`:
#'   - For sample size: `pA`, `pB`, `alpha`, `beta`
#'   - For power: `pA`, `pB`, `alpha`, `nB`
#'
#' - `"non-inferiority"`/`"equivalence"`:
#'   - For sample size: `pA`, `pB`, `delta`,  `alpha`, `beta`
#'   - For power: `pA`, `pB`, `delta`, `alpha`, `nB`
#'
#' @examples
#' # Sample size for `"equality"` test
#' or_size(pA = 0.4, pB = 0.25, kappa = 1,
#'         alpha = 0.05, beta = 0.2, test_type = "equality")
#'
#' # Power of `"equality"` test
#' or_size(pA = 0.4, pB = 0.25, kappa = 1,
#'         alpha = 0.05, nB = 156, test_type = "equality")
#'
#' # Sample size for `"non-inferiority"` test
#' or_size(pA = 0.4, pB = 0.25, delta = 0.2, kappa = 1,
#'         alpha = 0.05, beta = 0.2, test_type = "non-inferiority")
#'
#' # Power of `"non-inferiority"` test
#' or_size(pA = 0.4, pB = 0.25, delta = 0.2, kappa = 1,
#'         alpha = 0.05, nB = 242, test_type = "non-inferiority")
#'
#' # Sample size for `"equivalence"` test
#' or_size(pA = 0.25, pB = 0.25, delta = 0.5, kappa = 1,
#'         alpha = 0.05, beta = 0.2, test_type = "equivalence")
#'
#' # Power of `"equivalence"` test
#' or_size(pA = 0.25, pB = 0.25, delta = 0.5, kappa = 1,
#'         alpha = 0.05, nB = 366, test_type = "equivalence")
#'
#' @export
or_size <- function(pA, pB, delta = NULL, kappa = 1, alpha, beta = NULL, nB = NULL, test_type = "equality") {
  if (!is.null(beta)) {
    if (test_type == "equality") {
      return(ceiling((1/(kappa*pA*(1-pA))+1/(pB*(1-pB)))*((qnorm(1-alpha/2)+qnorm(1-beta))/log(pA*(1-pB)/pB/(1-pA)))^2))
    } else if (test_type == "non-inferiority") {
      return(ceiling((1/(kappa*pA*(1-pA))+1/(pB*(1-pB)))*((qnorm(1-alpha)+qnorm(1-beta))/(log(pA*(1-pB)/pB/(1-pA))-delta))^2))
    } else if (test_type == "equivalence") {
      return(ceiling((1/(kappa*pA*(1-pA))+1/(pB*(1-pB)))*((qnorm(1-alpha)+qnorm(1-beta/2))/(abs(log(pA*(1-pB)/pB/(1-pA)))-delta))^2))
    }
  } else if (!is.null(nB)) {
    if (test_type == "equality") {
      return(pnorm(log(pA*(1-pB)/pB/(1-pA))*sqrt(nB)/sqrt(1/(kappa*pA*(1-pA))+1/(pB*(1-pB)))-qnorm(1-alpha/2))+pnorm(-log(pA*(1-pB)/pB/(1-pA))*sqrt(nB)/sqrt(1/(kappa*pA*(1-pA))+1/(pB*(1-pB)))-qnorm(1-alpha/2)))
    } else if (test_type == "non-inferiority") {
      return(pnorm((log(pA*(1-pB)/pB/(1-pA))-delta)*sqrt(nB)/sqrt(1/(kappa*pA*(1-pA))+1/(pB*(1-pB)))-qnorm(1-alpha))+pnorm(-(log(pA*(1-pB)/pB/(1-pA))-delta)*sqrt(nB)/sqrt(1/(kappa*pA*(1-pA))+1/(pB*(1-pB)))-qnorm(1-alpha)))
    } else if (test_type == "equivalence") {
      return(2*(pnorm((abs(log(pA*(1-pB)/pB/(1-pA)))-delta)*sqrt(nB)/sqrt(1/(kappa*pA*(1-pA))+1/(pB*(1-pB)))-qnorm(1-alpha))+pnorm(-(abs(log(pA*(1-pB)/pB/(1-pA)))-delta)*sqrt(nB)/sqrt(1/(kappa*pA*(1-pA))+1/(pB*(1-pB)))-qnorm(1-alpha)))-1)
    }
  }
}
