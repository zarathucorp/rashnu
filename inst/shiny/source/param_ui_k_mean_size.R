param_ui_k_mean_size <- function() {
  tagList(
    selectInput("test_type", "Test Type", choices = c(
      "2-sided" = "2-side",
      "1-sided" = "1-side"
    )),
    radioButtons("mode", "Calculate", choices = c("Sample Size" = "size", "Power" = "power"), inline = TRUE),
    fluidRow(
      column(6, numericInput("muA", "muA", NULL)),
      column(6, numericInput("muB", "muB", NULL))
    ),
    conditionalPanel(
      "input.test_type == '1-side'",
      numericInput("kappa", "kappa", 1)
    ),
    conditionalPanel(
      "input.test_type == '2-side'",
      numericInput("sd", "sd", NULL)
    ),
    conditionalPanel(
      "input.test_type == '1-side'",
      fluidRow(
        column(6, numericInput("sdA", "sdA", NULL)),
        column(6, numericInput("sdB", "sdB", NULL))
      )
    ),
    numericInput("tau", "tau", NULL),
    fluidRow(
      column(6, numericInput("alpha", "alpha", 0.05)),
      conditionalPanel(
        "input.mode == 'size'",
        column(6, numericInput("beta", "beta", 0.2))
      ),
      conditionalPanel(
        "input.mode == 'power' && input.test_type == '2-side'",
        column(6, numericInput("n", "n", NULL))
      ),
      conditionalPanel(
        "input.mode == 'power' && input.test_type == '1-side'",
        column(6, numericInput("nA", "nA", NULL))
      )
    )
  )
}
