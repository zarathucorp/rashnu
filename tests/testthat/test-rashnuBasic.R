test_that("rashnuBasic returns a shiny.appobj", {
  app <- rashnuBasic()
  expect_s3_class(app, "shiny.appobj")
})
