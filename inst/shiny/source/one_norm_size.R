#' Sample Size or Power Calculation for One-Sample Normal Mean Test
#'
#' Calculates sample size or power for a two-sample normal mean test.
#'
#' @param mu Numeric. True mean.
#' @param mu0 Numeric. Null hypothesis mean.
#' @param sd Numeric. Standard deviation.
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
#' - For sample size: `"mu"`, `"mu0"`, `"sd"`, `"alpha"`, `"beta"`
#' - For power: `"mu"`, `"mu0"`, `"sd"`, `"alpha"`, `"n"`
#'
#' @examples
#' # Sample size
#' one_norm_size(mu = 2, mu0 = 1.5, sd = 1,
#'               alpha = 0.05, beta = 0.2)
#'
#' # Power
#' one_norm_size(mu = 2, mu0 = 1.5, sd = 1,
#'               alpha = 0.05, n = 32)
#'
#' @export
one_norm_size <- function(mu, mu0, sd, alpha, beta = NULL, n = NULL) {
  if (!is.null(beta)) {
    return(ceiling((sd*(qnorm(1-alpha/2)+qnorm(1-beta))/(mu-mu0))^2))
  } else if (!is.null(n)) {
    return(pnorm((mu-mu0)/sd*sqrt(n)-qnorm(1-alpha/2))+pnorm(-(mu-mu0)/sd*sqrt(n)-qnorm(1-alpha/2)))
  }
}
