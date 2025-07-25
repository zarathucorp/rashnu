main_ui_one_bino_size <- function() {
  withMathJax(
    tagList(
      tags$div(
        style = "font-size: 24px; font-weight: bold;",
        textOutput("result")
      ),
      hr(),
      h4("Calculate Sample Size Needed to Other: 1-Sample Binomial"),
      p("This calculator is useful for tests concerning whether a proportion, $p$, is equal to a reference value, $p_0$. The Null and Alternative hypotheses are"),
      p("$$H_0: p = p_0$$"),
      p("$$H_1: p \\ne p_0$$"),
      h4("Formulas"),
      p("This calculator uses the following formulas to compute sample size and power, respectively:"),
      p("$$n = p(1 - p) \\left( \\frac{z_{1 - \\alpha/2} + z_{1 - \\beta}}{p - p_0} \\right)^2$$"),
      p("$$1 - \\beta = \\Phi\\left( \\frac{p - p_0}{\\sqrt{\\frac{p(1 - p)}{n}}} - z_{1 - \\alpha/2} \\right) + \\Phi\\left( -\\frac{p - p_0}{\\sqrt{\\frac{p(1 - p)}{n}}} - z_{1 - \\alpha/2} \\right)$$"),
      p("where"),
      p("\\(n\\) is sample size"),
      p("\\(p_0\\) is the comparison value"),
      p("\\(\\Phi\\) is the standard Normal distribution function"),
      p("\\(\\Phi^{-1}\\) is the standard Normal quantile function"),
      p("\\(\\alpha\\) is Type I error"),
      p("\\(\\beta\\) is Type II error, meaning \\(1 - \\beta\\) is power")
    )
  )
}
