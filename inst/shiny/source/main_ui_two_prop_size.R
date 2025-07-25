main_ui_two_prop_size <- function(test_type) {
  withMathJax(
    tagList(
      tags$div(
        style = "font-size: 24px; font-weight: bold;",
        textOutput("result")
      ),
      hr(),
      switch(test_type,
             "2-side" = tagList(
               h4("Calculate Sample Size Needed to Compare 2 Proportions: 2-Sample, 2-Sided Equality"),
               p("This calculator is useful for tests concerning whether the proportions in two groups are different. Suppose the two groups are 'A' and 'B', and we collect a sample from both groups -- i.e. we have two samples. We perform a two-sample test to determine whether the proportion in group A, $p_A$, is different from the proportion in group B, $p_B$. The hypotheses are"),
               p("$$H_0:p_A-p_B=0$$"),
               p("$$H_1:p_A-p_B\\neq0$$"),
               p("where the ratio between the sample sizes of the two groups is"),
               p("$$\\kappa=\\frac{n_A}{n_B}$$"),
               h4("Formulas"),
               p("This calculator uses the following formulas to compute sample size and power, respectively:"),
               p("$$ n_A=\\kappa n_B \\;\\text{ and }\\; n_B=\\left(\\frac{p_A(1-p_A)}{\\kappa}+p_B(1-p_B)\\right) \\left(\\frac{z_{1-\\alpha/2}+z_{1-\\beta}}{p_A-p_B}\\right)^2$$"),
               p("$$1-\\beta= \\Phi\\left(z-z_{1-\\alpha/2}\\right)+\\Phi\\left(-z-z_{1-\\alpha/2}\\right) \\quad ,\\quad z=\\frac{p_A-p_B}{\\sqrt{\\frac{p_A(1-p_A)}{n_A}+\\frac{p_B(1-p_B)}{n_B}}}$$"),
               p("where"),
               p("\\(\\kappa=n_A/n_B\\) is the matching ratio"),
               p("\\(\\Phi\\) is the standard Normal distribution function"),
               p("\\(\\Phi^{-1}\\) is the standard Normal quantile function"),
               p("\\(\\alpha\\) is Type I error"),
               p("\\(\\beta\\) is Type II error, meaning \\(1 - \\beta\\) is power")
             ),
             "1-side" = tagList(
               h4("Calculate Sample Size Needed to Compare 2 Proportions: 2-Sample, 1-Sided"),
               p("This calculator is useful for tests concerning whether the proportions in two groups are different. Suppose the two groups are 'A' and 'B', and we collect a sample from both groups -- i.e. we have two samples. We perform a two-sample test to determine whether the proportion in group A, \\(p_A\\), is different from the proportion in group B, \\(p_B\\). The hypotheses are"),
               p("$$H_0:p_A=p_B$$"),
               p("$$H_1:p_A\\lt p_B$$"),
               p("or"),
               p("$$H_0:p_A=p_B$$"),
               p("$$H_1:p_A\\gt p_B$$"),
               p("where the ratio between the sample sizes of the two groups is"),
               p("$$\\kappa=\\frac{n_A}{n_B}$$"),
               h4("Formulas"),
               p("This calculator uses the following formulas to compute sample size and power, respectively:"),
               p("$$ n_A=\\kappa n_B \\;\\text{ and }\\; n_B=\\left(\\frac{p_A(1-p_A)}{\\kappa}+p_B(1-p_B)\\right) \\left(\\frac{z_{1-\\alpha}+z_{1-\\beta}}{p_A-p_B}\\right)^2$$"),
               p("$$1-\\beta=\\Phi\\left(\\frac{|p_A-p_B|}{\\sqrt{\\frac{p_A(1-p_A)}{n_A}+\\frac{p_B(1-p_B)}{n_B}}}-z_{1-\\alpha}\\right)$$"),
               p("where"),
               p("\\(\\kappa=n_A/n_B\\) is the matching ratio"),
               p("\\(\\Phi\\) is the standard Normal distribution function"),
               p("\\(\\Phi^{-1}\\) is the standard Normal quantile function"),
               p("\\(\\alpha\\) is Type I error"),
               p("\\(\\beta\\) is Type II error, meaning \\(1 - \\beta\\) is power")
             ),
             "non-inferiority" = tagList(
               h4("Calculate Sample Size Needed to Compare 2 Proportions: 2-Sample Non-Inferiority or Superiority"),
               p("This calculator is useful for the types of tests known as non-inferiority and superiority tests. Whether the null hypothesis represents 'non-inferiority' or 'superiority' depends on the context and whether the non-inferiority/superiority margin, \\(\\delta\\), is positive or negative. In this setting, we wish to test whether the proportion in group 'A', \\(p_A\\), is non-inferior/superior to the proportion in group 'B', \\(p_B\\). We collect a sample from both groups, and thus will conduct a two-sample test. The idea is that statistically significant differences between the proportions may not be of interest unless the difference is greater than a threshold, \\(\\delta\\). This is particularly popular in clinical studies, where the margin is chosen based on clinical judgement and subject-domain knowledge. The hypotheses to test are"),
               p("$$H_0:p_A-p_B\\le\\delta$$"),
               p("$$H_1:p_A-p_B>\\delta$$"),
               p("where \\(\\delta\\) is the superiority or non-inferiority margin and the ratio between the sample sizes of the two groups is"),
               p("$$\\kappa=\\frac{n_A}{n_B}$$"),
               h4("Formulas"),
               p("This calculator uses the following formulas to compute sample size and power, respectively:"),
               p("$$ n_A=\\kappa n_B \\;\\text{ and }\\; n_B=\\left(\\frac{p_A(1-p_A)}{\\kappa}+p_B(1-p_B)\\right) \\left(\\frac{z_{1-\\alpha}+z_{1-\\beta}}{p_A-p_B-\\delta}\\right)^2$$"),
               p("$$1-\\beta= \\Phi\\left(z-z_{1-\\alpha/2}\\right)+\\Phi\\left(-z-z_{1-\\alpha/2}\\right) \\quad ,\\quad z=\\frac{p_A-p_B-\\delta}{\\sqrt{\\frac{p_A(1-p_A)}{n_A}+\\frac{p_B(1-p_B)}{n_B}}}$$"),
               p("where"),
               p("\\(\\kappa=n_A/n_B\\) is the matching ratio"),
               p("\\(\\Phi\\) is the standard Normal distribution function"),
               p("\\(\\Phi^{-1}\\) is the standard Normal quantile function"),
               p("\\(\\alpha\\) is Type I error"),
               p("\\(\\beta\\) is Type II error, meaning \\(1 - \\beta\\) is power"),
               p("\\(\\delta\\) is the testing margin")
             ),
             "equivalence" = tagList(
               h4("Calculate Sample Size Needed to Compare 2 Proportions: 2-Sample Equivalence"),
               p("This calculator is useful when we wish to test whether the proportions in two groups are equivalent, without concern of which group's proportion is larger. Suppose we collect a sample from a group 'A' and a group 'B'; that is we collect two samples, and will conduct a two-sample test. For example, we may wish to test whether a new product is equivalent to an existing, industry standard product. Here, the 'burden of proof', so to speak, falls on the new product; that is, equivalence is actually represented by the alternative, rather than the null hypothesis."),
               p("$$H_0:|p_A-p_B|\\ge\\delta$$"),
               p("$$H_1:|p_A-p_B|<\\delta$$"),
               p("where \\(\\delta\\) is the superiority or non-inferiority margin and the ratio between the sample sizes of the two groups is"),
               p("$$\\kappa=\\frac{n_A}{n_B}$$"),
               h4("Formulas"),
               p("This calculator uses the following formulas to compute sample size and power, respectively:"),
               p("$$ n_A=\\kappa n_B \\;\\text{ and }\\; n_B=\\left(\\frac{p_A(1-p_A)}{\\kappa}+p_B(1-p_B)\\right) \\left(\\frac{z_{1-\\alpha}+z_{1-\\beta/2}}{|p_A-p_B|-\\delta}\\right)^2$$"),
               p("$$1-\\beta= 2\\left[\\Phi\\left(z-z_{1-\\alpha}\\right)+\\Phi\\left(-z-z_{1-\\alpha}\\right)\\right]-1 \\quad ,\\quad z=\\frac{|p_A-p_B|-\\delta}{\\sqrt{\\frac{p_A(1-p_A)}{n_A}+\\frac{p_B(1-p_B)}{n_B}}}$$"),
               p("where"),
               p("\\(\\kappa=n_A/n_B\\) is the matching ratio"),
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
