% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lakatosSampleSize.R
\name{lakatosSampleSize}
\alias{lakatosSampleSize}
\title{Sample Size Calculation Using Lakatos Method for Survival Analysis}
\usage{
lakatosSampleSize(
  syear,
  yrsurv1,
  yrsurv2,
  alloc,
  accrualTime,
  followTime,
  alpha,
  power,
  method = c("logrank", "gehan", "tarone-ware"),
  side = c("two.sided", "one.sided"),
  b = 24
)
}
\arguments{
\item{syear}{Survival time horizon in years.}

\item{yrsurv1}{Survival probability of the standard group at \code{syear}.}

\item{yrsurv2}{Survival probability of the test group at \code{syear}.}

\item{alloc}{Allocation ratio (Test / Standard). For equal allocation, use 1.}

\item{accrualTime}{Accrual period duration.}

\item{followTime}{Additional follow-up time after last patient is accrued.}

\item{alpha}{Significance level (e.g., 0.05 for two-sided tests).}

\item{power}{Desired statistical power (e.g., 0.8).}

\item{method}{Weighting method for test statistic. One of \code{"logrank"}, \code{"gehan"}, or \code{"tarone-ware"}.}

\item{side}{Type of test: \code{"two.sided"} or \code{"one.sided"}.}

\item{b}{Number of time divisions per year for numerical integration (default = 24).}
}
\value{
A list containing:
\describe{
\item{Sample_size_of_standard_group}{Required sample size in the standard group.}
\item{Sample_size_of_test_group}{Required sample size in the test group.}
\item{Total_sample_size}{Total sample size.}
\item{Expected_event_numbers_of_standard_group}{Expected number of events in the standard group.}
\item{Expected_event_numbers_of_test_group}{Expected number of events in the test group.}
\item{Total_expected_event_numbers}{Total number of expected events.}
\item{Actual_power}{Achieved power given the calculated sample size.}
\item{error}{(Optional) Error message when sample size cannot be calculated.}
}
}
\description{
Computes the required sample size and expected event numbers for two-group survival analysis
using Lakatos' method under exponential survival assumptions and varying weight functions (log-rank, Gehan, Tarone-Ware).
}
\examples{
lakatosSampleSize(
  syear = 2,
  yrsurv1 = 0.7,
  yrsurv2 = 0.6,
  alloc = 1,
  accrualTime = 1,
  followTime = 1,
  alpha = 0.05,
  power = 0.8,
  method = "logrank",
  side = "two.sided"
)


}
\references{
Lakatos E. (1988). Sample sizes based on the log-rank statistic in complex clinical trials.
\emph{Biometrics}, 44, 229–241.

Lakatos E, Lan KK. (1992). A comparison of sample size methods for the logrank statistic.
\emph{Statistics in Medicine}, 11(2), 179–191.

Web calculator (Superiority):
https://nshi.jp/en/js/twosurvyr/
}
