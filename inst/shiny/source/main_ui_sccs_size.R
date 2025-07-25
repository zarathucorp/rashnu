main_ui_sccs_size <- function() {
  withMathJax(
    tagList(
      tags$div(
        style = "font-size: 24px; font-weight: bold;",
        textOutput("result")
      ),
      hr(),
      h4("Calculate Sample Size Needed to Test Relative Incidence in Self Controlled Case Series Studies: SCCS, Alt-2"),
      p("The Self-Controlled Case Series (SCCS) method was originally developed by Farrington (1995) to compare relative incidence of adverse events following vaccination. In brief, the method compares incidence in a 'risk' time period shortly following exposure (vaccination) to the remainder of the observation period, the control period. Since it's initial development the SCCS method has been expanded in several ways, and has been used in a wide variety of pharmacoepidemiology studies. A salient feature of the method is that factors that do not vary with time (e.g. sex, race) are inherently accounted for."),
      p("Suppose each individual spends \\(r\\) proportion of the observation period in the exposed time period, and let \\(\\rho=e^{\\gamma}\\) represent the relative incidence; i.e. incidence is increased by a factor of \\(\\rho\\) in the exposed period compared to the control period. The hypotheses to be tested are"),
      p("$$H_0:\\rho=1$$"),
      p("$$H_1:\\rho \\neq 1$$"),
      p("This calculator implements the second method (i.e. alternative 2) for computing sample size for SCCS studies presented by Musonda et al. (2006). These computation formulas are based on the sampling distribution of \\(\\gamma=\\log(\\rho)\\)."),
      p("Under the null hypothesis \\(\\gamma\\) is approximately distributed"),
      p("$$ N\\left(0\\;,\\;\\frac{1}{n r (1-r)}\\right)$$"),
      p("where \\(n\\) is the number of individuals exposed during the observation period. Under the alternative hypothesis \\(\\gamma\\) is approximately distributed"),
      p("$$ N\\left(\\gamma\\;,\\;\\frac{(pr+1-r)^2}{npr(1-r)}\\right)$$"),
      p("This leads to the below formulas for power and sample size."),
      h4("Formulas"),
      p("This calculator uses the following formulas to compute sample size and power, respectively:"),
      p("$$n=\\left(\\frac{z_{1-\\alpha/2}+z_{1-\\beta}(\\rho r + 1 - r)/\\sqrt{\\rho}}{\\sqrt{r(1-r)}\\log(\\rho)}\\right)^2$$"),
      p("$$1-\\beta= \\Phi\\left( \\frac{\\gamma\\sqrt{n\\rho r(1-r)} - z_{1-\\alpha/2}\\sqrt{\\rho}}{\\rho r + 1 - r} \\right) $$"),
      p("where"),
      p("\\(n\\) is sample size"),
      p("\\(\\Phi\\) is the standard Normal distribution function"),
      p("\\(\\Phi^{-1}\\) is the standard Normal quantile function"),
      p("\\(\\alpha\\) is Type I error"),
      p("\\(\\beta\\) is Type II error, meaning \\(1 - \\beta\\) is power")
    )
  )
}
