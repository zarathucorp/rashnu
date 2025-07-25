main_ui_one_prop_size <- function(test_type) {
  withMathJax(
    tagList(
      tags$div(
        style = "font-size: 24px; font-weight: bold;",
        textOutput("result")
      ),
      hr(),
      switch(test_type,
             "2-side" = tagList(
               h4("Calculate Sample Size Needed to Test 1 Proportion: 1-Sample, 2-Sided Equality"),
               p("This calculator is useful for tests concerning whether a proportion, \\(p\\), is equal to a reference value, \\(p_0\\). The Null and Alternative hypotheses are"),
               p("$$H_0:p=p_0$$"),
               p("$$H_1:p\\neq p_0$$"),
               h4("Formulas"),
               p("This calculator uses the following formulas to compute sample size and power, respectively:"),
               p("$$n=p(1-p)\\left(\\frac{z_{1-\\alpha/2}+z_{1-\\beta}}{p-p_0}\\right)^2$$"),
               p("$$1-\\beta= \\Phi\\left(z-z_{1-\\alpha/2}\\right)+\\Phi\\left(-z-z_{1-\\alpha/2}\\right) \\quad ,\\quad z=\\frac{p-p_0}{\\sqrt{\\frac{p(1-p)}{n}}}$$"),
               p("where"),
               p("\\(n\\) is sample size"),
               p("\\(p_0\\) is the comparison value"),
               p("\\(\\Phi\\) is the standard Normal distribution function"),
               p("\\(\\Phi^{-1}\\) is the standard Normal quantile function"),
               p("\\(\\alpha\\) is Type I error"),
               p("\\(\\beta\\) is Type II error, meaning \\(1 - \\beta\\) is power")
             ),
             "1-side" = tagList(
               h4("Calculate Sample Size Needed to Test 1 Proportion: 1-Sample, 1-Sided"),
               p("This calculator is useful for tests concerning whether a proportion, \\(p\\), is equal to a reference value, \\(p_0\\). The Null and Alternative hypotheses are"),
               p("$$H_0:p=p_0$$"),
               p("$$H_1:p\\lt p_0$$"),
               p("or"),
               p("$$H_0:p=p_0$$"),
               p("$$H_1:p\\gt p_0$$"),
               h4("Formulas"),
               p("This calculator uses the following formulas to compute sample size and power, respectively:"),
               p("$$n=p_0(1-p_0)\\left(\\frac{z_{1-\\alpha}+z_{1-\\beta}\\sqrt{\\frac{p(1-p)}{p_0(1-p_0)}}}{p-p_0}\\right)^2$$"),
               p("$$1-\\beta=\\Phi\\left(\\sqrt{\\frac{p_0(1-p_0)}{p(1-p)}}\\left(\\frac{|p-p_0|\\sqrt{n}}{\\sqrt{p_0(1-p_0)}}-z_{1-\\alpha})\\right)\\right)$$"),
               p("where"),
               p("\\(n\\) is sample size"),
               p("\\(p_0\\) is the comparison value"),
               p("\\(\\Phi\\) is the standard Normal distribution function"),
               p("\\(\\Phi^{-1}\\) is the standard Normal quantile function"),
               p("\\(\\alpha\\) is Type I error"),
               p("\\(\\beta\\) is Type II error, meaning \\(1 - \\beta\\) is power")
             ),
             "non-inferiority" = tagList(
               h4("Calculate Sample Size Needed to Test 1 Proportion: 1-Sample Non-Inferiority or Superiority"),
               p("This calculator is useful for the types of tests known as non-inferiority and superiority tests. Whether the null hypothesis represents 'non-inferiority' or 'superiority' depends on the context and whether the non-inferiority/superiority margin, \\(\\delta\\), is positive or negative. In this setting, we wish to test whether a proportion, \\(p\\), is non-inferior/superior to a reference value, \\(p_0\\). The idea is that statistically significant differences between the proportion and the reference value may not be of interest unless the difference is greater than a threshold, \\(\\delta\\). This is particularly popular in clinical studies, where the margin is chosen based on clinical judgement and subject-domain knowledge. The hypotheses to test are"),
               p("$$H_0:p-p_0\\le\\delta$$"),
               p("$$H_1:p-p_0>\\delta$$"),
               p("and \\(\\delta\\) is the superiority or non-inferiority margin."),
               h4("Formulas"),
               p("This calculator uses the following formulas to compute sample size and power, respectively:"),
               p("$$n=p(1-p)\\left(\\frac{z_{1-\\alpha}+z_{1-\\beta}}{p-p_0-\\delta}\\right)^2$$"),
               p("$$1-\\beta= \\Phi\\left(z-z_{1-\\alpha}\\right)+\\Phi\\left(-z-z_{1-\\alpha}\\right) \\quad ,\\quad z=\\frac{p-p_0-\\delta}{\\sqrt{\\frac{p(1-p)}{n}}}$$"),
               p("where"),
               p("\\(n\\) is sample size"),
               p("\\(p_0\\) is the comparison value"),
               p("\\(\\Phi\\) is the standard Normal distribution function"),
               p("\\(\\Phi^{-1}\\) is the standard Normal quantile function"),
               p("\\(\\alpha\\) is Type I error"),
               p("\\(\\beta\\) is Type II error, meaning \\(1 - \\beta\\) is power"),
               p("\\(\\delta\\) is the testing margin")
             ),
             "equivalence" = tagList(
               h4("Calculate Sample Size Needed to Test 1 Proportion: 1-Sample Equivalence"),
               p("This calculator is useful when we wish to test whether a proportion, \\(p\\), is different from a gold standard reference value, \\(p_0\\). For example, we may wish to test whether a new product is equivalent to an existing, industry standard product. Here, the 'burden of proof', so to speak, falls on the new product; that is, equivalence is actually represented by the alternative, rather than the null hypothesis."),
               p("$$H_0:|p-p_0|\\ge\\delta$$"),
               p("$$H_1:|p-p_0|<\\delta$$"),
               h4("Formulas"),
               p("This calculator uses the following formulas to compute sample size and power, respectively:"),
               p("$$n=p(1-p)\\left(\\frac{z_{1-\\alpha}+z_{1-\\beta/2}}{|p-p_0|-\\delta}\\right)^2$$"),
               p("$$1-\\beta= 2\\left[\\Phi\\left(z-z_{1-\\alpha}\\right)+\\Phi\\left(-z-z_{1-\\alpha}\\right)\\right]-1 \\quad ,\\quad z=\\frac{|p-p_0|-\\delta}{\\sqrt{\\frac{p(1-p)}{n}}}$$"),
               p("where"),
               p("\\(n\\) is sample size"),
               p("\\(p_0\\) is the comparison value"),
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
