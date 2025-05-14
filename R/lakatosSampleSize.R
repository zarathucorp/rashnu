#' Sample Size Calculation Using Lakatos Method for Survival Analysis
#'
#' Computes the required sample size and expected event numbers for two-group survival analysis
#' using Lakatos' method under exponential survival assumptions and varying weight functions (log-rank, Gehan, Tarone-Ware).
#'
#' @param syear Survival time horizon in years.
#' @param yrsurv1 Survival probability of the standard group at `syear`.
#' @param yrsurv2 Survival probability of the test group at `syear`.
#' @param alloc Allocation ratio (Test / Standard). For equal allocation, use 1.
#' @param accrualTime Accrual period duration.
#' @param followTime Additional follow-up time after last patient is accrued.
#' @param alpha Significance level (e.g., 0.05 for two-sided tests).
#' @param power Desired statistical power (e.g., 0.8).
#' @param method Weighting method for test statistic. One of `"logrank"`, `"gehan"`, or `"tarone-ware"`.
#' @param side Type of test: `"two.sided"` or `"one.sided"`.
#' @param b Number of time divisions per year for numerical integration (default = 24).
#'
#' @return A list containing:
#' \describe{
#'   \item{Sample_size_of_standard_group}{Required sample size in the standard group.}
#'   \item{Sample_size_of_test_group}{Required sample size in the test group.}
#'   \item{Total_sample_size}{Total sample size.}
#'   \item{Expected_event_numbers_of_standard_group}{Expected number of events in the standard group.}
#'   \item{Expected_event_numbers_of_test_group}{Expected number of events in the test group.}
#'   \item{Total_expected_event_numbers}{Total number of expected events.}
#'   \item{Actual_power}{Achieved power given the calculated sample size.}
#'   \item{error}{(Optional) Error message when sample size cannot be calculated.}
#' }
#'
#' @examples
#' lakatosSampleSize(
#'   syear = 2,
#'   yrsurv1 = 0.7,
#'   yrsurv2 = 0.6,
#'   alloc = 1,
#'   accrualTime = 1,
#'   followTime = 1,
#'   alpha = 0.05,
#'   power = 0.8,
#'   method = "logrank",
#'   side = "two.sided"
#' )
#'
#' @import shiny
#' @import DT
#' @export
lakatosSampleSize <- function(syear, yrsurv1, yrsurv2, alloc,
                              accrualTime, followTime,
                              alpha, power,
                              method = c("logrank", "gehan", "tarone-ware"),
                              side = c("two.sided", "one.sided"),
                              b = 24
) {

  h1 <- -log(yrsurv1) / syear
  h2 <- -log(yrsurv2) / syear
  beta <- 1 - power


  method <- match.arg(method)
  side <- match.arg(side)
  totalTime <- accrualTime + followTime
  m <- floor(totalTime * b)

  ti <- seq(0, totalTime, length.out = m)
  n1 <- numeric(m)
  n2 <- numeric(m)

  allocRatio <- 1 / (1 + alloc)
  hr <- h2 / h1
  numer <- denom <- phi <- di <- wi <- e <- n <- 0
  eEvt1 <- eEvt2 <- 0

  #i in seq_along(ti)
  for (i in 1:m) {
    if (i == 1) {
      n1[i] <- allocRatio
      n2[i] <- 1 - allocRatio
    } else {
      if (ti[i - 1] <= followTime) {
        n1[i] <- n1[i - 1] * (1 - h1 / b)
        n2[i] <- n2[i - 1] * (1 - h2 / b)
      } else {
        n1[i] <- n1[i - 1] * (1 - h1 / b - 1 / (b * (totalTime - ti[i - 1])))
        n2[i] <- n2[i - 1] * (1 - h2 / b - 1 / (b * (totalTime - ti[i - 1])))
      }
    }

    phi <- n2[i] / n1[i]
    di <- (n1[i] * h1 + n2[i] * h2) / b

    # if(i<=10 | i>=1045) print(c(i, n1[i]+n2[i]))

    # if(1){
    #   print(i)
    #   print(n1[i]+n2[i])
    # }

    wi <- switch(method,
                 logrank = 1,
                 gehan = n1[i] + n2[i],
                 `tarone-ware` = sqrt(max( n1[i] + n2[i], 0)))

    numer <- numer + di * wi * ((hr * phi) / (1 + hr * phi) - phi / (1 + phi))
    denom <- denom + di * wi^2 * (phi / (1 + phi)^2)
    eEvt1 <- eEvt1 + n1[i] * h1 / b
    eEvt2 <- eEvt2 + n2[i] * h2 / b
  }

  e <- numer / sqrt(denom)
  result <- list()

  if (is.finite(e) && abs(e) > 0) {
    target_beta <- 1 - power

    repeat {
      n <- n + 1
      if (side == "two.sided") {
        pow <- pnorm(-sqrt(n) * e - qnorm(1 - alpha / 2)) +
          pnorm(sqrt(n) * e - qnorm(1 - alpha / 2))
      } else {
        pow <- pnorm(-sqrt(n) * e - qnorm(1 - alpha))
      }
      if (pow >= power) break
    }

    std_n <- ceiling(n * allocRatio)
    test_n <- ceiling(std_n * alloc)
    total_n <- std_n + test_n

    result <- list(
      Sample_size_of_standard_group = std_n,
      Sample_size_of_test_group = test_n,
      Total_sample_size = total_n,
      Expected_event_numbers_of_standard_group = round(n * eEvt1, 1),
      Expected_event_numbers_of_test_group = round(n * eEvt2, 1),
      Total_expected_event_numbers = round(n * (eEvt1 + eEvt2), 1),
      Actual_power = round(pow, 3)
    )
  } else {
    result <- list(
      error = "Non-finite or zero effect size detected. Unable to compute sample size."
    )
  }

  return(result)
}
