#' Sample Size or Power Calculation for K proportion
#'
#' Calculates sample size or power for a multiple-sample proportion test.
#'
#' @param pA Numeric. True proportion of group A.
#' @param pB Numeric. True proportion of group B.
#' @param tau Integer. Number of comparisons.
#' @param alpha Numeric. Type I error rate.
#' @param beta Numeric (optional). Type II error rate. Required for sample size calculation.
#' @param n Integer (optional). Sample size. Required when calculating power.
#'
#' @return Returns sample size (if `beta` is given), or power (if `n` is given).
#'
#' @note
#' Only one of `beta` (for sample size calculation) or `n` (for power calculation) should be specified.
#'
#' Required arguments:
#' - For sample size: `pA`, `pB`, `tau`, `alpha`, `beta`
#' - For power: `pA`, `pB`, `tau`, `alpha`, `n`
#'
#' @examples
#' # Sample size
#' k_prop_size(pA = 0.2, pB = 0.4, tau = 2,
#'             alpha = 0.05, beta = 0.2)
#'
#' # Power
#' k_prop_size(pA = 0.2, pB = 0.4, tau = 2,
#'             alpha = 0.05, n = 96)
#'
#' @export
k_prop_size <- function(pA, pB, tau, alpha, beta = NULL, n = NULL) {
  if (!is.null(beta)) {
    return(ceiling((pA*(1-pA)+pB*(1-pB))*((qnorm(1-alpha/2/tau)+qnorm(1-beta))/(pA-pB))^2))
  } else if (!is.null(n)) {
    return(pnorm((pA-pB)/sqrt(pA*(1-pA)/n+pB*(1-pB)/n)-qnorm(1-alpha/2/tau))+pnorm(-(pA-pB)/sqrt(pA*(1-pA)/n+pB*(1-pB)/n)-qnorm(1-alpha/2/tau)))
  }
}
