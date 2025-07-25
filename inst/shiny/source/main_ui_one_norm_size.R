main_ui_one_norm_size <- function() {
  withMathJax(
    tagList(
      tags$div(
        style = "font-size: 24px; font-weight: bold;",
        textOutput("result")
      ),
      hr(),
      h4("Calculate Sample Size Needed to Other: 1-Sample Normal"),
      p("This calculator is useful for tests concerning whether a mean, \\(\\mu\\), is equal to a reference value, \\(\\mu_0\\). The Null and Alternative hypotheses are"),
      p("$$H_0:\\mu=\\mu_0$$"),
      p("$$H_1:\\mu\\neq\\mu_0$$"),
      h4("Formulas"),
      p("This calculator uses the following formulas to compute sample size and power, respectively:"),
      p("$$n=\\left(\\sigma\\frac{z_{1-\\alpha/2}+z_{1-\\beta}}{\\mu-\\mu_0}\\right)^2$$"),
      p("$$1-\\beta=\\Phi\\left(\\frac{\\mu-\\mu_0}{\\sigma/\\sqrt{n}}-z_{1-\\alpha/2}\\right)+\\Phi\\left(-\\frac{\\mu-\\mu_0}{\\sigma/\\sqrt{n}}-z_{1-\\alpha/2}\\right)$$"),
      p("where"),
      p("where"),
      p("\\(n\\) is sample size"),
      p("\\(\\sigma\\) is standard deviation"),
      p("\\(\\Phi\\) is the standard Normal distribution function"),
      p("\\(\\Phi^{-1}\\) is the standard Normal quantile function"),
      p("\\(\\alpha\\) is Type I error"),
      p("\\(\\beta\\) is Type II error, meaning \\(1 - \\beta\\) is power")
    )
  )
}
