param_ui_or_size <- function() {
  tagList(
    selectInput("test_type", "Test Type", choices = c(
      "Equality" = "equality",
      "Non-inferiority" = "non-inferiority",
      "Equivalence" = "equivalence"
    )),
    radioButtons("mode", "Calculate", choices = c("Sample Size" = "size", "Power" = "power")),
    numericInput("pA", "pA", NULL),
    numericInput("pB", "pB", NULL),
    conditionalPanel(
      "input.test_type == 'non-inferiority' || input.test_type == 'equivalence'",
      numericInput("delta", "delta", NULL)
    ),
    numericInput("kappa", "Kappa", 1),
    numericInput("alpha", "alpha", 0.05),
    conditionalPanel(
      "input.mode == 'size'",
      numericInput("beta", "beta", 0.2)
    ),
    conditionalPanel(
      "input.mode == 'power'",
      numericInput("nB", "nB", NULL)
    )
  )
}
