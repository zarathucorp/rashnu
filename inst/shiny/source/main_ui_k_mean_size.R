main_ui_k_mean_size <- function(test_type) {
  withMathJax(
    tagList(
      tags$div(
        style = "font-size: 24px; font-weight: bold;",
        textOutput("result")
      ),
      hr(),
      switch(test_type,
             "2-side" = tagList(
               h4("Calculate Sample Size Needed to Compare k Means: 1-Way ANOVA Pairwise, 2-Sided Equality"),
               p("This calculator is useful for tests concerning whether the means of several groups are equal. The statistical model is called an Analysis of Variance, or ANOVA model. This calculator is for the particular situation where we wish to make", em("pairwise"), "comparisons between groups. That is, we test for equality between two groups at a time, and we make several of these comparisons."),
               p("For example, suppose we want to compare the means of three groups called foo, bar, and ack. These groups may represent groups of people that have been exposed to three different medical procedures, marketing schemes, etc. The complete list of pairwise comparisons are foo vs. bar, foo vs. ack, and bar vs. ack."),
               p("In more general terms, we may have \\(k\\) groups, meaning there are a total of \\(K \\equiv \\binom{k}{2} = k(k - 1)/2\\) possible pairwise comparisons. When we test \\(\\tau \\le K\\) of these pairwise comparisons, we have \\(\\tau\\) hypotheses of the form"),
               p("$$H_0:\\mu_A=\\mu_B$$"),
               p("$$H_1:\\mu_A\\ne\\mu_B$$"),
               p("where \\(\\mu_A\\) and \\(\\mu_B\\) represent the means of two of the \\(k\\) groups, groups 'A' and 'B'. We'll compute the required sample size for each of the \\(\\tau\\) comparisons, and total sample size needed is the largest of these. In the formula below, \\(n\\) represents the sample size in any one of these \\(\\tau\\) comparisons; that is, there are \\(n/2\\) people in the 'A' group, and \\(n/2\\) people in the 'B' group."),
               h4("Formulas"),
               p("This calculator uses the following formulas to compute sample size and power, respectively:"),
               p("$$n = 2\\left( \\sigma \\frac{z_{1 - \\alpha / (2\\tau)} + z_{1 - \\beta}}{\\mu_A - \\mu_B} \\right)^2$$"),
               p("$$1 - \\beta = \\Phi\\left( z - z_{1 - \\alpha / (2\\tau)} \\right) + \\Phi\\left( -z - z_{1 - \\alpha / (2\\tau)} \\right)\\quad ,\\quad z = \\frac{\\mu_A - \\mu_B}{\\sigma \\sqrt{\\frac{2}{n}}}$$"),
               p("where"),
               p("\\(n\\) is sample size"),
               p("\\(\\sigma\\) is standard deviation"),
               p("\\(\\Phi\\) is the standard Normal distribution function"),
               p("\\(\\Phi^{-1}\\) is the standard Normal quantile function"),
               p("\\(\\alpha\\) is Type I error"),
               p("\\(\\tau\\) is the number of comparisons to be made"),
               p("\\(\\beta\\) is Type II error, meaning \\(1 - \\beta\\) is power")
             ),
             "1-side" = tagList(
               h4("Calculate Sample Size Needed to Compare k Means: 1-Way ANOVA Pairwise, 1-Sided"),
               p("This calculator is useful for testing the means of several groups. The statistical model is called an Analysis of Variance, or ANOVA model. This calculator is for the particular situation where we wish to make pairwise comparisons between groups. That is, we compare two groups at a time, and we make several of these comparisons."),
               p("For example, suppose we want to compare the means of three groups called foo, bar, and ack. These groups may represent groups of people that have been exposed to three different medical procedures, marketing schemes, etc. The complete list of pairwise comparisons are foo vs. bar, foo vs. ack, and bar vs. ack."),
               p("In more general terms, we may have \\(k\\) groups, meaning there are a total of \\(K \\equiv \\binom{k}{2} = k(k - 1)/2\\) possible pairwise comparisons. When we test \\(\\tau \\le K\\) of these pairwise comparisons, we have \\(\\tau\\) hypotheses of the form"),
               p("$$H_0:\\mu_A = \\mu_B$$"),
               p("$$H_1:\\mu_A < \\mu_B$$"),
               p("or"),
               p("$$H_0:\\mu_A = \\mu_B$$"),
               p("$$H_1:\\mu_A < \\mu_B$$"),
               p("where \\(\\mu_A\\) and \\(\\mu_B\\) represent the means of two of the \\(k\\) groups, groups 'A' and 'B'. We'll compute the required sample size for each of the \\(\\tau\\) comparisons, and total sample size needed is the largest of these."),
               h4("Formulas"),
               p("This calculator uses the following formulas to compute sample size and power, respectively:"),
               p("$$n_A = \\left( \\sigma_A^2 + \\frac{\\sigma_B^2}{\\kappa} \\right) \\left( \\frac{z_{1 - \\alpha / \\tau} + z_{1 - \\beta}}{\\mu_A - \\mu_B} \\right)^2$$"),
               p("$$n_B = \\kappa\\, n_A$$"),
               p("$$1 - \\beta = \\Phi\\left( \\frac{|\\mu_A - \\mu_B| \\sqrt{n_A}}{\\sqrt{\\sigma_A^2 + \\sigma_B^2 / \\kappa}} - z_{1 - \\alpha / \\tau} \\right)$$"),
               p("where"),
               p("\\(\\kappa = n_A / n_B\\) is the matching ratio"),
               p("\\(\\sigma\\) is standard deviation"),
               p("\\(\\sigma_A\\) is standard deviation in Group \"A\""),
               p("\\(\\sigma_B\\) is standard deviation in Group \"B\""),
               p("\\(\\Phi\\) is the standard Normal distribution function"),
               p("\\(\\Phi^{-1}\\) is the standard Normal quantile function"),
               p("\\(\\alpha\\) is Type I error"),
               p("\\(\\tau\\) is the number of comparisons to be made"),
               p("\\(\\beta\\) is Type II error, meaning \\(1 - \\beta\\) is power")
             )
      )
    )
  )
}
