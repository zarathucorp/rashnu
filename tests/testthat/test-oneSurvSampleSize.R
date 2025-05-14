test_that("oneSurvSampleSize returns valid result with default log-log method", {
  res <- oneSurvSampleSize(
    survTime = 2,
    p1 = 0.75,
    p2 = 0.6,
    accrualTime = 1,
    followTime = 1,
    alpha = 0.05,
    power = 0.8,
    side = "two.sided",
    method = "log-log"
  )

  expect_type(res, "double")
  expect_named(res, c("SampleSize", "Power"))
  expect_true(res["SampleSize"] > 0)

})

methods <- c("arcsin", "log-log", "logit", "log", "log-swog", "identity")

for (m in methods) {
  test_that(paste("oneSurvSampleSize works with method:", m), {
    res <- oneSurvSampleSize(
      survTime = 12,
      p1 = 0.305,
      p2 = 0.435,
      accrualTime = 24,
      followTime = 24,
      alpha = 0.05,
      power = 0.8,
      side = "two.sided",
      method = m
    )

    expect_type(res, "double")
    expect_named(res, c("SampleSize", "Power"))
    expect_true(is.finite(res["SampleSize"]))
    expect_true(res["SampleSize"] > 0)

  })
}


test_that("oneSurvSampleSize handles no effect (p1 == p2)", {
  res <- oneSurvSampleSize(
    survTime = 2,
    p1 = 0.6,
    p2 = 0.6,
    accrualTime = 1,
    followTime = 1,
    alpha = 0.05,
    power = 0.8,
    side = "two.sided",
    method = "log-log"
  )

  expect_true(res["SampleSize"] > 10000)  # 아주 크거나 무한대
})
