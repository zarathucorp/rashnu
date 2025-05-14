test_that("lakatosSampleSize returns valid result with reasonable input", {
  res <- lakatosSampleSize(
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

  expect_type(res, "list")
  expect_true(is.null(res$error))  # 에러 메시지가 없어야 함
  expect_named(res, c(
    "Sample_size_of_standard_group",
    "Sample_size_of_test_group",
    "Total_sample_size",
    "Expected_event_numbers_of_standard_group",
    "Expected_event_numbers_of_test_group",
    "Total_expected_event_numbers",
    "Actual_power"
  ))
  expect_gt(res$Total_sample_size, 0)
  expect_lte(res$Actual_power, 1)
  expect_gte(res$Actual_power, 0.8)
})

test_that("lakatosSampleSize returns error when hazard ratios are equal", {
  res <- lakatosSampleSize(
    syear = 2,
    yrsurv1 = 0.7,
    yrsurv2 = 0.7,  # 동일한 생존율 → 효과 없음
    alloc = 1,
    accrualTime = 1,
    followTime = 1,
    alpha = 0.05,
    power = 0.8,
    method = "logrank",
    side = "two.sided"
  )

  expect_type(res, "list")
  expect_true(!is.null(res$error))
  expect_match(res$error, "effect size", ignore.case = TRUE)
})


test_that("lakatosSampleSize handles perfect survival input", {
  res <- lakatosSampleSize(
    syear = 2,
    yrsurv1 = 1,
    yrsurv2 = 0.6,
    alloc = 1,
    accrualTime = 1,
    followTime = 1,
    alpha = 0.05,
    power = 0.8,
    method = "logrank",
    side = "two.sided"
  )

  expect_type(res, "list")
  expect_true(!is.null(res$error))
})

