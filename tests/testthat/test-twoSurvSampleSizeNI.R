test_that("twoSurvSampleSizeNI returns valid output for standard inputs", {
  res <- twoSurvSampleSizeNI(
    syear = 2,
    yrsurv1 = 0.7,
    yrsurv2 = 0.65,
    alloc = 1,
    accrualTime = 1,
    followTime = 1,
    alpha = 0.025,
    power = 0.8,
    margin = 1.3
  )

  expect_type(res, "list")
  expect_named(res, c(
    "Sample_size_of_standard_group",
    "Sample_size_of_test_group",
    "Total_sample_size",
    "Expected_event_numbers_of_standard_group",
    "Expected_event_numbers_of_test_group",
    "Total_expected_event_numbers"
  ))

  expect_true(res$Sample_size_of_standard_group > 0)
  expect_true(res$Sample_size_of_test_group > 0)
  expect_true(res$Total_sample_size >= res$Sample_size_of_standard_group + res$Sample_size_of_test_group - 1)
  expect_true(res$Total_expected_event_numbers > 0)
})


test_that("twoSurvSampleSizeNI handles very small effect sizes", {
  res <- twoSurvSampleSizeNI(
    syear = 2,
    yrsurv1 = 0.7,
    yrsurv2 = 0.695,  # 거의 비슷한 생존율
    alloc = 1,
    accrualTime = 1,
    followTime = 1,
    alpha = 0.025,
    power = 0.8,
    margin = 1.3
  )

  expect_true(res$Total_sample_size > 1000)  # 효과 작으면 샘플 수 커져야 함
})

test_that("twoSurvSampleSizeNI handles large non-inferiority margin", {
  res <- twoSurvSampleSizeNI(
    syear = 2,
    yrsurv1 = 0.7,
    yrsurv2 = 0.65,
    alloc = 1,
    accrualTime = 1,
    followTime = 1,
    alpha = 0.025,
    power = 0.8,
    margin = 2.0
  )

  expect_true(res$Total_sample_size < 500)  # margin이 클수록 샘플 수는 작아짐
})

