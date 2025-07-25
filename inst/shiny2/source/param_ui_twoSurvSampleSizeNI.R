param_ui_twoSurvSampleSizeNI <- function() {
  tagList(
    numericInput("syear", "Survival Time", 12),
    fluidRow(column(6, numericInput("yrsurv1", "Survival Probability (Standard Group)", 0.3)),
             column(6, numericInput("yrsurv2", "Survival Probability (Test Group)", 0.5))),
    numericInput("alloc", "Allocation Ratio", 1),
    fluidRow(column(6, numericInput("accrualTime", "Accrual Time", 24)),
             column(6, numericInput("followTime", "FollowTime", 24))),
    fluidRow(column(6, numericInput("alpha", "Significance Level (alpha)", 0.025)),
             column(6, numericInput("power", HTML("Power (1- Beta)<br>&nbsp;"), 0.8))),
    numericInput("margin", "Non-inferiority Margin", 1.3),
    selectInput("side", "Hypothesis", choices = c(
      "one.sided" = "one.sided"
    ))
  )
}
