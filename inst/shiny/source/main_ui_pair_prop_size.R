main_ui_pair_prop_size <- function(test_type) {
  withMathJax(
    tagList(
      tags$div(
        style = "font-size: 24px; font-weight: bold;",
        textOutput("result")
      ),
      hr(),
      switch(test_type,
             "2-side" = tagList(
               h4("Calculate Sample Size Needed to Compare Paired Proportions: McNemar's Z-test, 2-Sided Equality"),
               p("This calculator is useful for tests comparing paired proportions. Suppose that our sample consists of pairs of subjects, and that each pair contains a subject from group 'A' and a subject from group 'B'. Further suppose that we wish to compare the probability that an event occurs in group 'A' to that in group 'B'. Example study designs include matched case-control studies and cross-over studies. Conceptually, the data can be listed as in the following table."),
               tags$table(
                 tags$tbody(
                   tags$tr(
                     tags$td(""),
                     tags$td(""),
                     tags$td(colspan = 2, "Group 'B'")
                   ),
                   tags$tr(
                     tags$td(""),
                     tags$td(""),
                     tags$td("Success"),
                     tags$td("Failure")
                   ),
                   tags$tr(
                     tags$td(rowspan = 2, "Group 'A'"),
                     tags$td("Success"),
                     tags$td(withMathJax("$$n_{11}$$")),
                     tags$td(withMathJax("$$n_{10}$$"))
                   ),
                   tags$tr(
                     tags$td("Failure"),
                     tags$td(withMathJax("$$n_{01}$$")),
                     tags$td(withMathJax("$$n_{00}$$"))
                   )
                 )
               ),
               p("Here, \\(n_{ij}\\) represents the number of pairs having $i$ successes in Group 'A' and $j$ successes in Group 'B'. The corresponding proportions are denoted \\(p_{ij}\\), with table"),
               tags$table(
                 tags$tbody(
                   tags$tr(
                     tags$td(""),
                     tags$td(""),
                     tags$td(colspan = 2, "Group 'B'")
                   ),
                   tags$tr(
                     tags$td(""),
                     tags$td(""),
                     tags$td("Success"),
                     tags$td("Failure")
                   ),
                   tags$tr(
                     tags$td(rowspan = 2, "Group 'A'"),
                     tags$td("Success"),
                     tags$td(withMathJax("$$p_{11}$$")),
                     tags$td(withMathJax("$$p_{10}$$"))
                   ),
                   tags$tr(
                     tags$td("Failure"),
                     tags$td(withMathJax("$$p_{01}$$")),
                     tags$td(withMathJax("$$p_{00}$$"))
                   )
                 )
               ),
               p("Interest is in comparing the following hypotheses:"),
               p("\\(H_0\\): Both groups have the same success probability", style = "text-align: center;"),
               p("\\(H_1\\): The success probability is not equal between the Groups", style = "text-align: center;"),
               p("Mathematically, this can be represented as"),
               p("$$H_0:p_{10}=p_{01}$$"),
               p("$$H_1:p_{10}\\neq p_{01}$$"),
               p("In the formulas below, we use the notation that"),
               p("$$p_{disc}=p_{10}+p_{01}$$"),
               p("and"),
               p("$$p_{diff}=p_{10}-p_{01}$$"),
               h4("Formulas"),
               p("This calculator uses the following formulas to compute sample size and power, respectively:"),
               p("$$ n=\\left(\\frac{z_{1-\\alpha/2}\\sqrt{p_{disc}}+z_{1-\\beta}\\sqrt{p_{disc}-p_{diff}^2}}{p_{diff}}\\right)^2$$"),
               p("$$1-\\beta=\\Phi\\left(\\frac{p_{diff}\\sqrt{n}-z_{1-\\alpha/2}\\sqrt{p_{disc}}}{\\sqrt{p_{disc}-p_{diff}^2}}\\right)$$"),
               p("\\(n\\) is sample size"),
               p("\\(\\Phi\\) is the standard Normal distribution function"),
               p("\\(\\Phi^{-1}\\) is the standard Normal quantile function"),
               p("\\(\\alpha\\) is Type I error"),
               p("\\(\\beta\\) is Type II error, meaning \\(1 - \\beta\\) is power")
             ),
             "1-side" = tagList(
               h4("Calculate Sample Size Needed to Compare Paired Proportions: McNemar's Z-test, 1-Sided"),
               p("This calculator is useful for tests comparing paired proportions. Suppose that our sample consists of pairs of subjects, and that each pair contains a subject from group 'A' and a subject from group 'B'. Further suppose that we wish to compare the probability that an event occurs in group 'A' to that in group 'B'. Example study designs include matched case-control studies and cross-over studies. Conceptually, the data can be listed as in the following table."),
               tags$table(
                 tags$tbody(
                   tags$tr(
                     tags$td(""),
                     tags$td(""),
                     tags$td(colspan = 2, "Group 'B'")
                   ),
                   tags$tr(
                     tags$td(""),
                     tags$td(""),
                     tags$td("Success"),
                     tags$td("Failure")
                   ),
                   tags$tr(
                     tags$td(rowspan = 2, "Group 'A'"),
                     tags$td("Success"),
                     tags$td(withMathJax("$$n_{11}$$")),
                     tags$td(withMathJax("$$n_{10}$$"))
                   ),
                   tags$tr(
                     tags$td("Failure"),
                     tags$td(withMathJax("$$n_{01}$$")),
                     tags$td(withMathJax("$$n_{00}$$"))
                   )
                 )
               ),
               p("Here, \\(n_{ij}\\) represents the number of pairs having \\(i\\) successes in Group 'A' and \\(j\\) successes in Group 'B'. The corresponding proportions are denoted \\(p_{ij}\\), with table"),
               tags$table(
                 tags$tbody(
                   tags$tr(
                     tags$td(""),
                     tags$td(""),
                     tags$td(colspan = 2, "Group 'B'")
                   ),
                   tags$tr(
                     tags$td(""),
                     tags$td(""),
                     tags$td("Success"),
                     tags$td("Failure")
                   ),
                   tags$tr(
                     tags$td(rowspan = 2, "Group 'A'"),
                     tags$td("Success"),
                     tags$td(withMathJax("$$p_{11}$$")),
                     tags$td(withMathJax("$$n_{10}$$"))
                   ),
                   tags$tr(
                     tags$td("Failure"),
                     tags$td(withMathJax("$$p_{01}$$")),
                     tags$td(withMathJax("$$p_{00}$$"))
                   )
                 )
               ),
               p("Interest is in comparing the following hypotheses:"),
               p("\\(H_0\\): Both groups have the same success probability", style = "text-align: center;"),
               p("\\(H_1\\): The success probability is not equal between the Groups", style = "text-align: center;"),
               p("Mathematically, this can be represented as"),
               p("$$H_0:p_{10}=p_{01}$$"),
               p("$$H_1:p_{10}\\lt p_{01}$$"),
               p("or"),
               p("$$H_0:p_{10}=p_{01}$$"),
               p("$$H_1:p_{10}\\gt p_{01}$$"),
               p("In the formulas below, we use the notation that"),
               p("or"),
               p("$$p_{disc}=p_{10}+p_{01}$$"),
               p("and"),
               p("$$p_{diff}=p_{10}-p_{01}$$"),
               h4("Formulas"),
               p("This calculator uses the following formulas to compute sample size and power, respectively:"),
               p("$$ n=\\left(\\frac{z_{1-\\alpha}\\sqrt{p_{disc}}+z_{1-\\beta}\\sqrt{p_{disc}-p_{diff}^2}}{p_{diff}}\\right)^2$$"),
               p("$$1-\\beta=\\Phi\\left(\\frac{|p_{diff}|\\sqrt{n}-z_{1-\\alpha}\\sqrt{p_{disc}}}{\\sqrt{p_{disc}-p_{diff}^2}}\\right)$$"),
               p("where"),
               p("\\(n\\) is sample size"),
               p("\\(\\Phi\\) is the standard Normal distribution function"),
               p("\\(\\Phi^{-1}\\) is the standard Normal quantile function"),
               p("\\(\\alpha\\) is Type I error"),
               p("\\(\\beta\\) is Type II error, meaning \\(1 - \\beta\\) is power")
             )
      )
    )
  )
}
