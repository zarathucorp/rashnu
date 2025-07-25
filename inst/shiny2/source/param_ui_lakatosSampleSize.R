param_ui_lakatosSampleSize <- function() {
  tagList(
    numericInput("syear", "Survival Time", 12),
    fluidRow(column(6, numericInput("yearsurv1", "Survivial Probability (Standard Group)", 0.3)),
             column(6, numericInput("yearsurv2", "Survivial Probability (Test Group)", 0.5))),
    numericInput("alloc", "Allocation Ratio", 1),
    fluidRow(column(6, numericInput("accrualTime", "Accrual Time", 24)),
             column(6, numericInput("followTime", "Follow-up Time", 24))),
    fluidRow(column(6, numericInput("alpha", "Significance Level (alpha)", 0.025)),
             column(6, numericInput("power", HTML("Power (1- Beta)<br>&nbsp;"), 0.8))),
    selectInput("method", "Test Method", choices = c(
      "logrank" = "logrank",
      "gehan" = "gehan",
      "tarone-ware" = "tarone-ware"
    )),
    selectInput("side", "Hypothesis", choices = c(
      "two.sided" = "two.sided",
      "one.sided" = "one.sided"
    ))
  )
}
