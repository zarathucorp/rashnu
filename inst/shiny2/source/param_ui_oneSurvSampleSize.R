param_ui_oneSurvSampleSize <- function() {
  tagList(
    numericInput("survTime", "Survival Time", 12),
    fluidRow(column(6, numericInput("p2", "Null Survivial Probability", 0.3)),
             column(6, numericInput("p1", "Alternative Survivial Probability", 0.5))),
    fluidRow(column(6, numericInput("accrualTime", "Accrual Time", 24)),
             column(6, numericInput("followTime", "Follow-up Time", 24))),
    fluidRow(column(6, numericInput("alpha", "Significance Level (alpha)", 0.025)),
             column(6, numericInput("power", HTML("Power (1- Beta)<br>&nbsp;"), 0.8))),
    selectInput("side", "Hypothesis", choices = c(
      "two.sided" = "two.sided",
      "one.sided" = "one.sided"
    )),
    selectInput("method", "Transformation", choices = c(
      "arcsin" = "arcsin",
      "log-log" = "log-log",
      "logit" = "logit",
      "log" = "log",
      "log-swog" = "log-swog",
      "identity" = "identity"
    ))
  )
}
