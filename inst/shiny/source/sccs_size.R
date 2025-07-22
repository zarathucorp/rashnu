#' Sample Size or Power for Self-Controlled Case Series (SCCS)
#'
#' Calculates sample size or power for self-controlled case series studies.
#'
#' @param p Numeric. True relative incidence (risk period vs baseline).
#' @param r Numeric. Proportion of observation time that is risk period.
#' @param alpha Numeric. Type I error rate.
#' @param beta Numeric (optional). Type II error rate. Required for sample size calculation.
#' @param n Integer (optional). Sample size. Required for power calculation.
#'
#' @return Numeric. Returns sample size (if `beta` is given), or power (if `n` is given).
#'
#' @note
#' Only one of `beta` (for sample size calculation) or `n` (for power calculation) should be specified.
#'
#' Required arguments:
#' - For sample size: `p`, `r`, `alpha`, `beta`
#' - For power: `p`, `r`, `alpha`, `n`
#'
#' @examples
#' # Sample size
#' sccs_size(p = 3, r = 42/365,
#'           alpha = 0.05, beta = 0.2)
#'
#' # Power
#' sccs_size(p = 3, r = 42/365,
#'           alpha = 0.05, n = 54)
#'
#' @export
sccs_size <- function(p, r, alpha, beta = NULL, n = NULL) {
  if (!is.null(beta)) {
    return(ceiling((qnorm(1-alpha/2)+qnorm(1-beta)*(p*r+1-r)/sqrt(p))^2/(r*(1-r)*log(p)^2)))
  } else if (!is.null(n)) {
    return(pnorm((log(p)*sqrt(n*p*r*(1-r))-qnorm(1-alpha/2)*sqrt(p))/(p*r+1-r)))
  }
}
