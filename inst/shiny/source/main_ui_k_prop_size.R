main_ui_k_prop_size <- function() {
  withMathJax(
    tagList(
      tags$div(
        style = "font-size: 24px; font-weight: bold;",
        textOutput("result")
      ),
      hr(),
      h4("Calculate Sample Size Needed to Compare k Proportions: 1-Way ANOVA Pairwise"),
      p("This calculator is useful for tests concerning whether the proportions in several groups are equal. The statistical model is called an Analysis of Variance, or ANOVA model. This calculator is for the particular situation where we wish to make pairwise comparisons between groups. That is, we test for equality between two groups at a time, and we make several of these comparisons."),
      p("For example, suppose we want to compare the proportions in three groups called foo, bar, and ack. These groups may represent groups of people that have been exposed to three different medical procedures, marketing schemes, etc. The complete list of pairwise comparisons are foo vs. bar, foo vs. ack, and bar vs. ack."),
      p("In more general terms, we may have \\(k\\) groups, meaning there are a total of \\(K \\equiv \\binom{k}{2} = k(k - 1)/2\\) possible pairwise comparisons. When we test \\(\\tau \\le K\\) of these pairwise comparisons, we have \\(\\tau\\) hypotheses of the form"),
      p("$$H_0: p_A = p_B$$"),
      p("$$H_1: p_A \\ne p_B$$"),
      p("where \\(p_A\\) and \\(p_B\\) represent the proportions in two of the \\(k\\) groups, groups 'A' and 'B'. We'll compute the required sample size for each of the \\(\\tau\\) comparisons, and total sample size needed is the largest of these. In the formula below, \\(n\\) represents the sample size in any one of these \\(\\tau\\) comparisons; that is, there are \\(n/2\\) people in the 'A' group, and \\(n/2\\) people in the 'B' group."),
      h4("Formulas"),
      p("This calculator uses the following formulas to compute sample size and power, respectively:"),
      p("$$n = \\left( p_A(1 - p_A) + p_B(1 - p_B) \\right) \\left( \\frac{z_{1 - \\alpha / (2\\tau)} + z_{1 - \\beta}}{p_A - p_B} \\right)^2$$"),
      p("$$1 - \\beta = \\Phi\\left( z - z_{1 - \\alpha / (2\\tau)} \\right) + \\Phi\\left( -z - z_{1 - \\alpha / (2\\tau)} \\right)\\quad ,\\quad z = \\frac{p_A - p_B}{\\sqrt{\\frac{p_A(1 - p_A)}{n} + \\frac{p_B(1 - p_B)}{n}}}$$"),
      p("where"),
      p("\\(n\\) is sample size"),
      p("\\(\\Phi\\) is the standard Normal distribution function"),
      p("\\(\\Phi^{-1}\\) is the standard Normal quantile function"),
      p("\\(\\alpha\\) is Type I error"),
      p("\\(\\tau\\) is the number of comparisons to be made"),
      p("\\(\\beta\\) is Type II error, meaning \\(1 - \\beta\\) is power")
    )
  )
}
