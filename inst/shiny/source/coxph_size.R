#' Sample Size or Power for Cox Proportional Hazards
#'
#' Calculates sample size or power for a cox proportional hazards model.
#'
#' @param hr Numeric. True hazard ratio.
#' @param hr0 Numeric (optional). Null hypothesis hazard ratio. Required for `"2-side"`, `"non-inferiority"` test.
#' @param delta Numeric (optional). Margin for `"equivalence test"`. Required for `"equivalence"` test.
#' @param pE Numeric. Overall event probability.
#' @param pA Numeric. Proportion of group A.
#' @param alpha Numeric. Type I error rate.
#' @param beta Numeric (optional). Type II error rate. Required for sample size calculation.
#' @param n Integer (optional). Sample size. Required for power calculation.
#' @param test_type Character. `"2-side"`, `"non-inferiority"`, or `"equivalence"`. Default is `"2-side"`.
#'
#' @return Numeric. Returns sample size (if `beta` is given), or power (if `n` is given).
#'
#' @note
#' Only one of `beta` (for sample size calculation) or `n` (for power calculation) should be specified.
#'
#' Required arguments by `test_type`:
#' - `"2-side"`/`"non-inferiority`:
#'   - For sample size: `hr`, `hr0`, `pE`, `pA`, `alpha`, `beta`
#'   - For power: `hr`, `hr0`, `pE`, `pA`, `alpha`, `n`
#'
#' - `"equivalence"`:
#'   - For sample size: `hr`, `delta`, `pE`, `pA`, `alpha`, `beta`
#'   - For power: `hr`, `delta`, `pE`, `pA`, `alpha`, `n`
#'
#' @examples
#' # Sample size for a `"2-side"` test
#' coxph_size(hr = 2, hr0 = 1, pE = 0.8, pA = 0.5,
#'            alpha = 0.05, beta = 0.2, test_type = "2-side")
#'
#' # Power of `"2-side"` test
#' coxph_size(hr = 2, hr0 = 1, pE = 0.8, pA = 0.5,
#'            alpha = 0.05, n = 82, test_type = "2-side")
#'
#' # Sample size for `"non-inferiority"` test
#' coxph_size(hr = 2, hr0 = 1, pE = 0.8, pA = 0.5,
#'            alpha = 0.025, beta = 0.2, test_type = "non-inferiority")
#'
#' # Power of `"non-inferiority"` test
#' coxph_size(hr = 2, hr0 = 1, pE = 0.8, pA = 0.5,
#'            alpha = 0.025, n = 82, test_type = "non-inferiority")
#'
#' # Sample size for `"equivalence"` test
#' coxph_size(hr = 1, delta = 0.5, pE = 0.8, pA = 0.5,
#'            alpha = 0.05, beta = 0.2, test_type = "equivalence")
#'
#' # Power of `"equivalence"` test
#' coxph_size(hr = 1, delta = 0.5, pE = 0.8, pA = 0.5,
#'            alpha = 0.05, n = 172, test_type = "equivalence")
#'
#' @export
coxph_size <- function(hr, hr0 = NULL, delta = NULL, pE, pA, alpha, beta = NULL, n = NULL, test_type = "2-side") {
  if (!is.null(beta)) {
    if (test_type == "2-side") {
      return(ceiling(((qnorm(1-alpha/2)+qnorm(1-beta))/(log(hr)-log(hr0)))^2/(pA*(1-pA)*pE)))
    } else if (test_type == "non-inferiority") {
      return(ceiling(((qnorm(1-alpha)+qnorm(1-beta))/(log(hr)-log(hr0)))^2/(pA*(1-pA)*pE)))
    } else if (test_type == "equivalence") {
      return(ceiling(((qnorm(1-alpha)+qnorm(1-beta/2))/(delta-abs(log(hr))))^2/(pA*(1-pA)*pE)))
    }
  } else if (!is.null(n)) {
    if (test_type == "2-side") {
      return(pnorm((log(hr)-log(hr0))*sqrt(n*pA*(1-pA)*pE)-qnorm(1-alpha/2)))
    } else if (test_type == "non-inferiority") {
      return(pnorm((log(hr)-log(hr0))*sqrt(n*pA*(1-pA)*pE)-qnorm(1-alpha)))
    } else if (test_type == "equivalence") {
      return(2*pnorm((delta-abs(log(hr)))*sqrt(n*pA*(1-pA)*pE)-qnorm(1-alpha))-1)
    }
  }
}
