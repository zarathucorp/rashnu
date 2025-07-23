param_ui_coxph_size <- function() {
  tagList(
    selectInput("test_type", "Test Type", choices = c(
      "2-sided" = "2-side",
      "Non-inferiority" = "non-inferiority",
      "Equivalence" = "equivalence"
    )),
    radioButtons("mode", "Calculate", choices = c("Sample Size" = "size", "Power" = "power")),
    numericInput("hr", "hr", NULL),
    conditionalPanel(
      "input.test_type == '2-side' || input.test_type == 'non-inferiority'",
      numericInput("hr0", "hr0", NULL)
    ),
    conditionalPanel(
      "input.test_type == 'equivalence'",
      numericInput("delta", "delta", NULL)
    ),
    numericInput("pE", "pE", NULL),
    numericInput("pA", "pA", NULL),
    numericInput("alpha", "alpha", 0.05),
    conditionalPanel(
      "input.mode == 'size'",
      numericInput("beta", "beta", 0.2)
    ),
    conditionalPanel(
      "input.mode == 'power'",
      numericInput("n", "n", NULL)
    )
  )
}
