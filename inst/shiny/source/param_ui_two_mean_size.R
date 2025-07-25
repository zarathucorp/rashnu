param_ui_two_mean_size <- function() {
  tagList(
    selectInput("test_type", "Test Type", choices = c(
      "2-sided" = "2-side",
      "1-sided" = "1-side",
      "Non-inferiority" = "non-inferiority",
      "Equivalence" = "equivalence"
    )),
    radioButtons("mode", "Calculate", choices = c("Sample Size" = "size", "Power" = "power"), inline = TRUE),
    fluidRow(
      column(6, numericInput("muA", "muA", NULL)),
      column(6, numericInput("muB", "muB", NULL))
    ),
    conditionalPanel(
      "input.test_type == 'non-inferiority' || input.test_type == 'equivalence'",
      numericInput("delta", "delta", NULL)
    ),
    numericInput("kappa", "Kappa", 1),
    conditionalPanel(
      "input.test_type != '1-side'",
      numericInput("sd", "sd", NULL)
    ),
    conditionalPanel(
      "input.test_type == '1-side'",
      fluidRow(
        column(6, numericInput("sdA", "sdA", NULL)),
        column(6, numericInput("sdB", "sdB", NULL))
      )
    ),
    fluidRow(
      column(6, numericInput("alpha", "alpha", 0.05)),
      conditionalPanel(
        "input.mode == 'size'",
        column(6, numericInput("beta", "beta", 0.2))
      ),
      conditionalPanel(
        "input.mode == 'power' && input.test_type == '2-side'",
        column(6, numericInput("nB", "nB", NULL))
      ),
      conditionalPanel(
        "input.mode == 'power' && input.test_type == '1-side'",
        column(6, numericInput("nA", "nA", NULL))
      )
    )
  )
}
