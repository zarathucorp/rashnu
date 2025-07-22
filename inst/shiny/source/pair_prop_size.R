#' Sample Size or Power for Paired-Sample Proportion Test
#'
#' Calculates sample size or power for a paired-sample proportion test.
#'
#' @param p01 Numeric. Proportion of discordant pairs with (before = 1, after = 0).
#' @param p10 Numeric. Proportion of discordant pairs with (before = 0, after = 1).
#' @param alpha Numeric. Type I error rate.
#' @param beta Numeric (optional). Type II error rate. Required for sample size calculation.
#' @param n Integer (optional). Sample size. Required for power calculation.
#' @param test_type Character. `"2-side"` or `"1-side"`. Default is `"2-side"`.
#'
#' @return Numeric. Returns sample size (if `beta` is given), or power (if `n` is given).
#'
#' @note
#' Only one of `beta` (for sample size calculation) or `n` (for power calculation) should be specified.
#'
#' Required arguments:
#' - For sample size: `p01`, `p10`, `alpha`, `beta`
#' - For power: `p01`, `p10`, `alpha`, `n`
#'
#' @examples
#' # Sample size for `"2-side"` test
#' pair_prop_size(p01 = 0.45, p10 = 0.05,
#'                alpha = 0.1, beta = 0.1, test_type = "2-side")
#'
#' # Power of `"2-side"` test
#' pair_prop_size(p01 = 0.45, p10 = 0.05,
#'                alpha = 0.1, n = 23, test_type = "2-side")
#'
#' # Sample size for `"1-side"` test
#' pair_prop_size(p01 = 0.45, p10 = 0.05,
#'                alpha = 0.05, beta = 0.1, test_type = "1-side")
#'
#' # Power of `"1-side"` test
#' pair_prop_size(p01 = 0.45, p10 = 0.05,
#'                alpha = 0.05, n = 23, test_type = "1-side")
#'
#' @export
pair_prop_size <- function(p01, p10, alpha, beta = NULL, n = NULL, test_type = "2-side") {
  if (!is.null(beta)) {
    if (test_type == "2-side") {
      return(ceiling(((qnorm(1-alpha/2)*sqrt(p10+p01)+qnorm(1-beta)*sqrt(p10+p01-(p10-p01)^2))/(p10-p01))^2))
    } else if (test_type == "1-side") {
      return(ceiling(((qnorm(1-alpha)*sqrt(p10+p01)+qnorm(1-beta)*sqrt(p10+p01-(p10-p01)^2))/(p10-p01))^2))
    }
  } else if (!is.null(n)) {
    if (test_type == "2-side") {
      return(pnorm(((p10-p01)*sqrt(n)-qnorm(1-alpha/2)*sqrt(p10+p01))/sqrt(p10+p01-(p10-p01)^2))+pnorm((-(p10-p01)*sqrt(n)-qnorm(1-alpha/2)*sqrt(p10+p01))/sqrt(p10+p01-(p10-p01)^2)))
    } else if (test_type == "1-side") {
      return(pnorm((abs(p10-p01)*sqrt(n)-qnorm(1-alpha)*sqrt(p10+p01))/sqrt(p10+p01-(p10-p01)^2)))
    }
  }
}
