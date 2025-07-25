param_ui_pair_prop_size <- function() {
  tagList(
    selectInput("test_type", "Test Type", choices = c(
      "2-sided" = "2-side",
      "1-sided" = "1-side"
    )),
    radioButtons("mode", "Calculate", choices = c("Sample Size" = "size", "Power" = "power"), inline = TRUE),
    fluidRow(
      column(6, numericInput("p01", "p01", NULL)),
      column(6, numericInput("p10", "p10", NULL))
    ),
    fluidRow(
      column(6, numericInput("alpha", "alpha", 0.05)),
      column(6,
             conditionalPanel("input.mode == 'size'", numericInput("beta", "beta", 0.2)),
             conditionalPanel("input.mode == 'power'", numericInput("n", "n", NULL))
      )
    )
  )
}
