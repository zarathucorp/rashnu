#' Sample Size Calculation for Two-Group Non-Inferiority Survival Study
#'
#' Calculates the required sample size and expected event numbers for a non-inferiority trial with two survival curves,
#' using piecewise integration of hazard functions under exponential survival assumptions.
#'
#' @param syear Survival time horizon (e.g., median survival time) in years.
#' @param yrsurv1 Survival probability of the standard group at `syear`.
#' @param yrsurv2 Survival probability of the test group at `syear`.
#' @param alloc Allocation ratio (Test / Standard), e.g., 1 means equal allocation.
#' @param accrualTime Duration of patient accrual period.
#' @param followTime Follow-up period after last patient is accrued.
#' @param alpha One-sided significance level (e.g., 0.025).
#' @param power Desired statistical power (e.g., 0.8).
#' @param margin Non-inferiority margin for hazard ratio (HR).
#'
#' @return A list containing:
#' \describe{
#'   \item{Sample_size_of_standard_group}{Required sample size in the standard group.}
#'   \item{Sample_size_of_test_group}{Required sample size in the test group.}
#'   \item{Total_sample_size}{Total sample size.}
#'   \item{Expected_event_numbers_of_standard_group}{Expected number of events in the standard group.}
#'   \item{Expected_event_numbers_of_test_group}{Expected number of events in the test group.}
#'   \item{Total_expected_event_numbers}{Total number of expected events across both groups.}
#' }
#'
#' @examples
#' twoSurvSampleSizeNI(
#'   syear = 2,
#'   yrsurv1 = 0.7,
#'   yrsurv2 = 0.65,
#'   alloc = 1,
#'   accrualTime = 1,
#'   followTime = 1,
#'   alpha = 0.025,
#'   power = 0.8,
#'   margin = 1.3
#' )
#' @importFrom stats pnorm qnorm
#' @export
twoSurvSampleSizeNI <- function(syear, yrsurv1, yrsurv2, alloc, accrualTime, followTime,  alpha, power, margin) {

  h1 <- -log(yrsurv1) / syear
  h2 <- -log(yrsurv2) / syear
  beta <- 1 - power

  totalTime <- accrualTime + followTime
  hr1 <- h2 / h1
  hr0 <- margin

  p2 <- alloc / (1 + alloc)
  p1 <- 1 - p2

  za <- qnorm(1 - alpha)
  zb <- qnorm(1 - beta)

  nk <- 5000
  w <- totalTime / nk

  i0 <- (p1 * h1 + p2 * h2) / (p1 + hr0 * p2)^2 * 0.5 * w
  i1 <- (p1 * h1 + p2 * h2) / (p1 + hr1 * p2)^2 * 0.5 * w
  om <- (p1 * h1 + p2 * h2) / ((p1 + hr0 * p2) * (p1 + hr1 * p2)) * 0.5 * w
  d1 <- 0
  d2 <- 0

  for (k in 1:(nk - 1)) {
    t <- k * w
    s1 <- exp(-h1 * t)
    s2 <- s1^hr1
    f1 <- h1 * s1
    f2 <- hr1 * h1 * s2

    if (t <= followTime) {
      g <- 1.0
      dg <- 0.0
    } else {
      g <- -t / accrualTime + totalTime / accrualTime
      dg <- -1.0 / accrualTime
    }

    i0 <- i0 + g * s1 * s2 * (p1 * f1 + p2 * f2) / (p1 * s1 + hr0 * p2 * s2)^2 * w
    i1 <- i1 + g * s1 * s2 * (p1 * f1 + p2 * f2) / (p1 * s1 + hr1 * p2 * s2)^2 * w
    om <- om + g * s1 * s2 * (p1 * f1 + p2 * f2) / ((p1 * s1 + hr0 * p2 * s2) * (p1 * s1 + hr1 * p2 * s2)) * w
    d1 <- d1 + dg * s1 * w
    d2 <- d2 + dg * s2 * w
  }

  i0 <- hr0 * p1 * p2 * i0
  i1 <- hr1 * p1 * p2 * i1
  om <- (hr0 - hr1) * p1 * p2 * om

  n <- ((sqrt(i0) * za + sqrt(i1) * zb) / om)^2
  sample_std <- ceiling(n * p1)
  sample_test <- ceiling(sample_std * alloc)
  list(
    Sample_size_of_standard_group = sample_std,
    Sample_size_of_test_group = sample_test,
    Total_sample_size = sample_std + sample_test,
    Expected_event_numbers_of_standard_group = round(n * p1 * (1 + d1), 1),
    Expected_event_numbers_of_test_group = round(n * p2 * (1 + d2), 1),
    Total_expected_event_numbers = round(n * p1 * (1 + d1) + n * p2 * (1 + d2), 1)
  )

}
