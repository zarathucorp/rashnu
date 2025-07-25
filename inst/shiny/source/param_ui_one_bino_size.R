param_ui_one_bino_size <- function() {
  tagList(
    radioButtons("mode", "Calculate", choices = c("Sample Size" = "size", "Power" = "power"), inline = TRUE),
    fluidRow(
      column(6, numericInput("p", "p", NULL)),
      column(6, numericInput("p0", "p0", NULL))
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
