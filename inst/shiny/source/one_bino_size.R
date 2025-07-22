#' Sample Size or Power for One-Sample Binomial Proportion Test
#'
#' Calculates sample size or power for a two-sample binomial proportion test.
#'
#' @param p Numeric. True proportion.
#' @param p0 Numeric. Null hypothesis proportion.
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
#' - For sample size: `"p"`, `"p0"`, `"alpha"`, `"beta"`
#' - For power: `"p"`, `"p0"`, `"alpha"`, `"n"`
#'
#' @examples
#' # Required sample size
#' one_bino_size(p = 0.5, p0 = 0.3,
#'               alpha = 0.05, beta = 0.2)
#'
#' # Power
#' one_bino_size(p = 0.5, p0 = 0.3,
#'               alpha = 0.05, n = 50)
#'
#' @export
one_bino_size <- function(p, p0, alpha, beta = NULL, n = NULL) {
  if (!is.null(beta)) {
    return(ceiling(p*(1-p)*((qnorm(1-alpha/2)+qnorm(1-beta))/(p-p0))^2))
  } else if (!is.null(n)) {
    return(pnorm((p-p0)/sqrt(p*(1-p)/n)-qnorm(1-alpha/2))+pnorm(-(p-p0)/sqrt(p*(1-p)/n)-qnorm(1-alpha/2)))
  }
}
