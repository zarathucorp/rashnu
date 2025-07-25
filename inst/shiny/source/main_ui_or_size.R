main_ui_or_size <- function(test_type) {
  withMathJax(
    tagList(
      tags$div(
        style = "font-size: 24px; font-weight: bold;",
        textOutput("result")
      ),
      hr(),
      switch(test_type,
             "equality" = tagList(
               h4("Calculate Sample Size Needed to Test Odds Ratio: Equality"),
               p("This calculator is useful for tests concerning whether the odds ratio, \\(OR\\), between two groups is different from the null value of 1. Suppose the two groups are 'A' and 'B', and we collect a sample from both groups -- i.e. we have two samples. We perform a two-sample test to determine whether the odds of the outcome in group A, \\(p_A(1 - p_A)\\), is different from the odds of the outcome in group B, \\(p_B(1 - p_B)\\), where \\(p_A\\) and \\(p_B\\) are the probabilities of the outcome in the two groups. The hypotheses are"),
               p("$$H_0:OR=1$$"),
               p("$$H_1:OR\\neq1$$"),
               p("where the ratio between the sample sizes of the two groups is"),
               p("$$\\kappa=\\frac{n_A}{n_B}$$"),
               p("and"),
               p("$$OR=\\frac{p_A(1-p_B)}{p_B(1-p_A)}$$"),
               h4("Formulas"),
               p("This calculator uses the following formulas to compute sample size and power, respectively:"),
               p("$$ n_A=\\kappa n_B \\;\\text{ and }\\; n_B=\\left(\\frac{1}{\\kappa p_A(1-p_A)}+\\frac{1}{p_B(1-p_B)}\\right) \\left(\\frac{z_{1-\\alpha/2}+z_{1-\\beta}}{\\ln(OR)}\\right)^2$$"),
               p("$$1 - \\beta = \\Phi\\left( z - z_{1 - \\alpha/2} \\right) + \\Phi\\left( -z - z_{1 - \\alpha/2} \\right)\\quad ,\\quad z = \\frac{\\ln(OR) \\sqrt{n_B}}{\\sqrt{\\left( \\frac{1}{\\kappa p_A (1 - p_A)} + \\frac{1}{p_B (1 - p_B)} \\right)}}$$"),
               p("where"),
               p("$$OR=\\frac{p_A(1-p_B)}{p_B(1-p_A)}$$"),
               p("and where"),
               p("\\(\\kappa = n_A / n_B\\) is the matching ratio"),
               p("\\(\\Phi\\) is the standard Normal distribution function"),
               p("\\(\\Phi^{-1}\\) is the standard Normal quantile function"),
               p("\\(\\alpha\\) is Type I error"),
               p("\\(\\beta\\) is Type II error, meaning \\(1 - \\beta\\) is power")
             ),
             "non-inferiority" = tagList(
               h4("Calculate Sample Size Needed to Test Odds Ratio: Non-Inferiority or Superiority"),
               p("This calculator is useful for the types of tests known as non-inferiority and superiority tests. Whether the null hypothesis represents 'non-inferiority' or 'superiority' depends on the context and whether the non-inferiority/superiority margin, \\(\\delta\\), is positive or negative. In this setting, we wish to test whether the odds of an outcome in group 'A', \\(p_A(1-p_A)\\), is non-inferior/superior to the odds of the outcome in group 'B', \\(p_B(1-p_B)\\), where \\(p_A\\) and \\(p_B\\) are the probabilities of the outcome in the two groups. We collect a sample from both groups, and thus will conduct a two-sample test. The idea is that statistically significant differences between the proportions may not be of interest unless the difference is greater than a threshold. This is particularly popular in clinical studies, where the margin is chosen based on clinical judgement and subject-domain knowledge. The hypotheses to test are"),
               p("$$H_0:\\ln(OR)\\le\\delta$$"),
               p("$$H_1:\\ln(OR)>\\delta$$"),
               p("where \\(\\delta\\) is the superiority or non-inferiority margin on the log scale, and the ratio between the sample sizes of the two groups is"),
               p("$$\\kappa=\\frac{n_A}{n_B}$$"),
               h4("Formulas"),
               p("This calculator uses the following formulas to compute sample size and power, respectively:"),
               p("$$ n_A=\\kappa n_B \\;\\text{ and }\\; n_B=\\left(\\frac{1}{\\kappa p_A(1-p_A)}+\\frac{1}{p_B(1-p_B)}\\right) \\left(\\frac{z_{1-\\alpha}+z_{1-\\beta}}{\\ln(OR)-\\delta}\\right)^2$$"),
               p("$$1-\\beta= \\Phi\\left(z-z_{1-\\alpha}\\right)+\\Phi\\left(-z-z_{1-\\alpha}\\right) \\quad ,\\quad z=\\frac{(\\ln(OR)-\\delta)\\sqrt{n_B}}{\\sqrt{\\frac{1}{\\kappa p_A(1-p_A)}+\\frac{1}{p_B(1-p_B)}}}$$"),
               p("where"),
               p("$$OR=\\frac{p_A(1-p_B)}{p_B(1-p_A)}$$"),
               p("and where"),
               p("\\(\\kappa = n_A / n_B\\) is the matching ratio"),
               p("\\(\\Phi\\) is the standard Normal distribution function"),
               p("\\(\\Phi^{-1}\\) is the standard Normal quantile function"),
               p("\\(\\alpha\\) is Type I error"),
               p("\\(\\beta\\) is Type II error, meaning \\(1 - \\beta\\) is power"),
               p("\\(\\delta\\) is the testing margin")
             ),
             "equivalence" = tagList(
               h4("Calculate Sample Size Needed to Test Odds Ratio: Equivalence"),
               p("This calculator is useful when we wish to test whether the odds of an outcome in two groups are equivalent, without concern of which group's odds is larger. Suppose we collect a sample from a group 'A' and a group 'B'; that is we collect two samples, and will conduct a two-sample test. For example, we may wish to test whether a new product is equivalent to an existing, industry standard product. Here, the 'burden of proof', so to speak, falls on the new product; that is, equivalence is actually represented by the alternative, rather than the null hypothesis."),
               p("$$H_0:|\\ln(OR)|\\ge\\delta$$"),
               p("$$H_1:|\\ln(OR)|<\\delta$$"),
               p("where \\(\\delta\\) is the superiority or non-inferiority margin and the ratio between the sample sizes of the two groups is"),
               p("$$\\kappa=\\frac{n_A}{n_B}$$"),
               h4("Formulas"),
               p("This calculator uses the following formulas to compute sample size and power, respectively:"),
               p("$$ n_A=\\kappa n_B \\;\\text{ and }\\; n_B=\\left(\\frac{1}{\\kappa p_A(1-p_A)}+\\frac{1}{p_B(1-p_B)}\\right) \\left(\\frac{z_{1-\\alpha}+z_{1-\\beta/2}}{|\\ln(OR)|-\\delta}\\right)^2$$"),
               p("$$1-\\beta= \\Phi\\left(z-z_{1-\\alpha}\\right)+\\Phi\\left(-z-z_{1-\\alpha}\\right) \\quad ,\\quad z=\\frac{(|\\ln(OR)|-\\delta)\\sqrt{n_B}}{\\sqrt{\\frac{1}{\\kappa p_A(1-p_A)}+\\frac{1}{p_B(1-p_B)}}}$$"),
               p("where"),
               p("$$OR=\\frac{p_A(1-p_B)}{p_B(1-p_A)}$$"),
               p("and where"),
               p("\\(\\kappa = n_A / n_B\\) is the matching ratio"),
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
