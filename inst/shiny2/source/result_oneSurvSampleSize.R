result_oneSurvSampleSize <- function(input) {
  oneSurvSampleSize(
    survTime = input$survTime,
    p2 = input$p2,
    p1 = input$p1,
    accrualTime = input$accrualTime,
    followTime = input$followTime,
    alpha = input$alpha,
    power = input$power,
    side = input$side,
    method = input$method
  )
}
