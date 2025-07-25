result_twoSurvSampleSizeNI <- function(input) {
    twoSurvSampleSizeNI(
      syear = input$syear,
      yrsurv1 = input$yrsurv1,
      yrsurv2 = input$yrsurv2,
      alloc = input$alloc,
      accrualTime = input$accrualTime,
      followTime = input$followTime,
      alpha = input$alpha,
      power = input$power,
      margin = input$margin
    )
}
