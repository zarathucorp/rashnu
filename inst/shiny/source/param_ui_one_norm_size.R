param_ui_one_norm_size <- function() {
  tagList(
    radioButtons("mode", "Calculate", choices = c("Sample Size" = "size", "Power" = "power"), inline = TRUE),
    fluidRow(
      column(6, numericInput("mu", "mu", NULL)),
      column(6, numericInput("mu0", "mu0", NULL))
    ),
    numericInput("sd", "sd", NULL),
    fluidRow(
      column(6, numericInput("alpha", "alpha", 0.05)),
      column(6,
             conditionalPanel("input.mode == 'size'", numericInput("beta", "beta", 0.2)),
             conditionalPanel("input.mode == 'power'", numericInput("n", "n", NULL))
      )
    )
  )
}
