main_ui_one_mean_size <- function(test_type) {
  withMathJax(
    tagList(
      tags$div(
        style = "font-size: 24px; font-weight: bold;",
        textOutput("result")
      ),
      hr(),
      switch(test_type,
             "2-side" = tagList(
               h4("Calculate Sample Size Needed to Test 1 Mean: 1-Sample, 2-Sided Equality"),
               p("This calculator is useful for tests concerning whether a mean, \\(\\mu\\), is equal to a reference value, \\(\\mu_0\\). The Null and Alternative hypotheses are"),
               p("$$H_0:\\mu=\\mu_0$$"),
               p("$$H_1:\\mu\\neq\\mu_0$$"),
               h4("Formulas"),
               p("This calculator uses the following formulas to compute sample size and power, respectively:"),
               p("$$n = \\left( \\sigma \\frac{z_{1 - \\alpha/2} + z_{1 - \\beta}}{\\mu - \\mu_0} \\right)^2$$"),
               p("$$1 - \\beta = \\Phi\\left( z - z_{1 - \\alpha/2} \\right) + \\Phi\\left( -z - z_{1 - \\alpha/2} \\right)\\quad ,\\quad z = \\frac{\\mu - \\mu_0}{\\sigma / \\sqrt{n}}$$"),
               p("where"),
               p("\\(n\\) is sample size"),
               p("\\(\\sigma\\) is standard deviation"),
               p("\\(\\Phi\\) is the standard Normal distribution function"),
               p("\\(\\Phi^{-1}\\) is the standard Normal quantile function"),
               p("\\(\\alpha\\) is Type I error"),
               p("\\(\\beta\\) is Type II error, meaning \\(1 - \\beta\\) is power")
             ),
             "1-side" = tagList(
               h4("Calculate Sample Size Needed to Test 1 Mean: 1-Sample, 1-Sided"),
               p("This calculator is useful for tests concerning whether a mean, \\(\\mu\\), is equal to a reference value, \\(\\mu_0\\). The Null and Alternative hypotheses is either"),
               p("$$H_0:\\mu=\\mu_0$$"),
               p("$$H_1:\\mu\\lt\\mu_0$$"),
               p("or"),
               p("$$H_0:\\mu=\\mu_0$$"),
               p("$$H_1:\\mu\\gt\\mu_0$$"),
               h4("Formulas"),
               p("This calculator uses the following formulas to compute sample size and power, respectively:"),
               p("$$n=\\left(\\sigma\\frac{z_{1-\\alpha}+z_{1-\\beta}}{\\mu-\\mu_0}\\right)^2$$"),
               p("$$1-\\beta=\\Phi\\left(\\frac{|\\mu-\\mu_0|}{\\sigma/\\sqrt{n}}-z_{1-\\alpha}\\right)$$"),
               p("where"),
               p("\\(n\\) is sample size"),
               p("\\(\\sigma\\) is standard deviation"),
               p("\\(\\Phi\\) is the standard Normal distribution function"),
               p("\\(\\Phi^{-1}\\) is the standard Normal quantile function"),
               p("\\(\\alpha\\) is Type I error"),
               p("\\(\\beta\\) is Type II error, meaning \\(1 - \\beta\\) is power")
             ),
             "non-inferiority" = tagList(
               h4("Calculate Sample Size Needed to Test 1 Mean: 1-Sample Non-Inferiority or Superiority"),
               p("This calculator is useful for the types of tests known as non-inferiority and superiority tests. Whether the null hypothesis represents 'non-inferiority' or 'superiority' depends on the context and whether the non-inferiority/superiority margin, \\(\\delta\\), is positive or negative. In this setting, we wish to test whether a mean, \\(\\mu\\), is non-inferior/superior to a reference value, \\(\\mu_0\\). The idea is that statistically significant differences between the mean and the reference value may not be of interest unless the difference is greater than a threshold, \\(\\delta\\). This is particularly popular in clinical studies, where the margin is chosen based on clinical judgement and subject-domain knowledge. The hypotheses to test are"),
               p("$$H_0:\\mu-\\mu_0\\le\\delta$$"),
               p("$$H_1:\\mu-\\mu_0>\\delta$$"),
               p("and \\(\\delta\\) is the superiority or non-inferiority margin."),
               h4("Formulas"),
               p("This calculator uses the following formulas to compute sample size and power, respectively:"),
               p("$$n=\\left(\\sigma\\frac{z_{1-\\alpha}+z_{1-\\beta}}{\\mu-\\mu_0-\\delta}\\right)^2$$"),
               p("$$1-\\beta= \\Phi\\left(z-z_{1-\\alpha}\\right)+\\Phi\\left(-z-z_{1-\\alpha}\\right) \\quad ,\\quad z=\\frac{\\mu-\\mu_0-\\delta}{\\sigma/\\sqrt{n}}$$"),
               p("where"),
               p("\\(n\\) is sample size"),
               p("\\(\\sigma\\) is standard deviation"),
               p("\\(\\Phi\\) is the standard Normal distribution function"),
               p("\\(\\Phi^{-1}\\) is the standard Normal quantile function"),
               p("\\(\\alpha\\) is Type I error"),
               p("\\(\\beta\\) is Type II error, meaning \\(1 - \\beta\\) is power"),
               p("\\(\\delta\\) is the testing margin")
             ),
             "equivalence" = tagList(
               h4("Calculate Sample Size Needed to Test 1 Mean: 1-Sample Equivalence"),
               p("This calculator is useful when we wish to test whether a mean, \\(\\mu\\), is different from a gold standard reference value, \\(\\mu_0\\). For example, we may wish to test whether a new product is equivalent to an existing, industry standard product. Here, the 'burden of proof', so to speak, falls on the new product; that is, equivalence is actually represented by the alternative, rather than the null hypothesis."),
               p("$$H_0:|\\mu-\\mu_0|\\ge\\delta$$"),
               p("$$H_1:|\\mu-\\mu_0|<\\delta$$"),
               h4("Formulas"),
               p("This calculator uses the following formulas to compute sample size and power, respectively:"),
               p("$$n=\\left(\\sigma\\frac{z_{1-\\alpha}+z_{1-\\beta/2}}{\\delta-|\\mu-\\mu_0|}\\right)^2$$"),
               p("$$1-\\beta=2\\left[\\Phi\\left(z-z_{1-\\alpha}\\right)+\\Phi\\left(-z-z_{1-\\alpha}\\right)\\right]-1 \\quad ,\\quad z=\\frac{|\\mu-\\mu_0|-\\delta}{\\sigma/\\sqrt{n}}$$"),
               p("where"),
               p("\\(n\\) is sample size"),
               p("\\(\\sigma\\) is standard deviation"),
               p("\\(\\Phi\\) is the standard Normal distribution function"),
               p("\\(\\Phi^{-1}\\) is the standard Normal quantile function"),
               p("\\(\\alpha\\) is Type I error"),
               p("\\(\\beta\\) is Type II error, meaning \\(1 - \\beta\\) is power"),
               p("\\(\\delta\\) is the testing margin")
             )
      )
    )
  )
}
