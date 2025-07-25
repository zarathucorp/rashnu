#' One-Sample Survival Study Sample Size or Power Calculation
#'
#' Calculates the required sample size or power for a single-arm survival study using various transformation-based methods,
#' including arcsine-square root, log-log, logit, and others. This function assumes an exponential survival model.
#'
#' @param survTime Time point at which survival is evaluated (e.g., median follow-up time).
#' @param p1 Expected survival probability under the alternative hypothesis.
#' @param p2 Survival probability under the null hypothesis.
#' @param accrualTime Patient accrual period.
#' @param followTime Additional follow-up period after accrual ends.
#' @param alpha Significance level (e.g., 0.05).
#' @param power Desired statistical power (e.g., 0.8).
#' @param side Type of hypothesis test. Either `"two.sided"` (default) or `"one.sided"`.
#' @param method Transformation method for comparison. One of `"arcsin"`, `"log-log"`, `"logit"`, `"log"`, `"log-swog"`, `"identity"`.
#'
#' @return A named numeric vector with:
#' \describe{
#'   \item{SampleSize}{Calculated required sample size.}
#'   \item{Power}{Achieved power with the calculated sample size.}
#' }
#'
#' @examples
#' oneSurvSampleSize(
#'   survTime = 2,
#'   p1 = 0.75,
#'   p2 = 0.6,
#'   accrualTime = 1,
#'   followTime = 1,
#'   alpha = 0.05,
#'   power = 0.8,
#'   side = "two.sided",
#'   method = "log-log"
#' )
#'
#'
#' @references
#' Fleming TR, Harrington DP. (1991). *Counting Processes and Survival Analysis*.
#' New York: Wiley, pp. 236–237, Example 6.3.1.
#'
#' Andersen PK, Borgan O, Gill RD, Keiding N. (1993). *Statistical Models Based on Counting Processes*.
#' New York: Springer-Verlag, pp. 176–287, Section IV.1–3.
#'
#' Bie O, Borgan O, Liestol K. (1987). Confidence intervals and confidence bands for the cumulative hazard rate function and their small sample properties.
#' *Scandinavian Journal of Statistics*, 14(3), 221–233.
#'
#' Borgan O, Liestol K. (1990). A note on confidence intervals and bands for the survival function based on transformations.
#' *Scandinavian Journal of Statistics*, 17(1), 35–41.
#'
#' Nagashima K, Noma H, Sato Y, Gosho M. (2020). Sample size calculations for single-arm survival studies
#' using transformations of the Kaplan–Meier estimator. *Pharmaceutical Statistics*. https://doi.org/10.1002/pst.2090
#' Available at: https://arxiv.org/abs/2012.03355
#'
#' Web calculator (One-sample):
#' https://nshi.jp/en/js/onesurvyr/
#'
#' @importFrom stats pnorm qnorm
#' @export
oneSurvSampleSize <- function(
    survTime, p1, p2,accrualTime, followTime,
    alpha, power, side = c("two.sided", "one.sided"), method = c("arcsin", "log-log", "logit", "log", "log-swog", "identity")
){
  h1 = -log(p1) / survTime
  h2 = -log(p2) / survTime

  beta <- 1 - power


  a <- accrualTime
  f <- followTime
  b <- a + f
  s <- survTime
  nk <- 5000
  result <- numeric(2)

  if (side == "one.sided") {
    za <- qnorm(1 - alpha)
  } else {
    za <- qnorm(1 - alpha / 2)
  }
  zb <- qnorm(1 - beta)

  if (s < f) {
    avar1 <- exp(h1 * s) - 1
    avar2 <- exp(h2 * s) - 1
  } else {
    w <- (s - f) / nk
    t_seq <- seq(f + w, s - w, length.out = nk - 1)
    avar1 <- exp(h1 * f) / (b - f) * 0.5 + sum(exp(h1 * t_seq) / (b - t_seq)) + exp(h1 * s) / (b - s) * 0.5
    avar2 <- exp(h2 * f) / (b - f) * 0.5 + sum(exp(h2 * t_seq) / (b - t_seq)) + exp(h2 * s) / (b - s) * 0.5

    avar1 <- w * avar1 * h1 * a + exp(h1 * f) - 1
    avar2 <- w * avar2 * h2 * a + exp(h2 * f) - 1
  }

  method <- tolower(method)
  if (method == "arcsin") {
    ncp <- abs(asin(sqrt(p1)) - asin(sqrt(p2)))
    avar2 <- avar2 * exp(-h2 * s) / (1 - exp(-h2 * s)) * 0.25
    n <- ceiling(avar2 * ((za + zb) / ncp)^2)
    power <- pnorm(-za + ncp * sqrt(n) / sqrt(avar2))
  } else if (method == "log-log") {
    ncp <- abs(log(-log(p1)) - log(-log(p2)))
    avar2 <- avar2 / (h2^2 * s^2)
    n <- ceiling(avar2 * ((za + zb) / ncp)^2)
    power <- pnorm(-za + ncp * sqrt(n) / sqrt(avar2))
  } else if (method == "logit") {
    ncp <- abs(log(p1 / (1 - p1)) - log(p2 / (1 - p2)))
    avar2 <- avar2 / (1 - exp(-h2 * s))^2
    n <- ceiling(avar2 * ((za + zb) / ncp)^2)
    power <- pnorm(-za + ncp * sqrt(n) / sqrt(avar2))
  } else if (method == "log") {
    ncp <- abs(log(p1) - log(p2))
    n <- ceiling(avar2 * ((za + zb) / ncp)^2)
    power <- pnorm(-za + ncp * sqrt(n) / sqrt(avar2))
  } else if (method == "log-swog") {
    ncp <- abs(log(p1) - log(p2))
    n <- ceiling(((sqrt(avar2) * za + sqrt(avar1) * zb) / ncp)^2)
    power <- pnorm(-za * sqrt(avar2) / sqrt(avar1) + ncp * sqrt(n) / sqrt(avar1))
  } else {
    # identity
    ncp <- abs(p1 - p2)
    avar2 <- avar2 * exp(-2 * h2 * s)
    n <- ceiling(avar2 * ((za + zb) / ncp)^2)
    power <- pnorm(-za + ncp * sqrt(n) / sqrt(avar2))
  }

  result[1] <- n
  result[2] <- round(power,3)
  names(result) <- c("SampleSize", "Power")
  return(result)
}
