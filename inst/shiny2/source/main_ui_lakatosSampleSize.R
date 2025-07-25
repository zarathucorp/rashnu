main_ui_lakatosSampleSize <- function() {
  withMathJax(
    tagList(
      tags$div(
        style = "font-size: 22px; font-weight: bold;",
        tableOutput("result")
      ),
      hr(),
      h4("References"),
      p("Jung SH, Chow SC. On sample size calculation for comparing survival curves under general hypothesis testing. Journal of Biopharmaceutical Statistics 2012; 22(3):485-495."),
      p("Lakatos E. Sample sizes based on the log-rank statistic in complex clinical trials. Biometrics 1988; 44:229-241."),
      p("Lakatos E, Lan KK. A comparison of sample size methods for the logrank statistic. Statistics in Medicine 1992; 11(2):179-191."),
      p("Fleming TR, Harrington DP. Counting Processes and Survival Analysis. New York: Wiley, 1991, 236-237, Example 6.3.1."),
      p("Andersen PK, Borgan Ø, Gill RD, Keiding N. Statistical Models Based on Counting Processes. New York: Springer-Verlag, 1993, 176-287, Section IV.1-3."),
      p("Bie O, Borgan Ø, LiestØl K. Confidence intervals and confidence bands for the cumulative hazard rate function and their small sample properties. Scandinavian Journal of Statistics 1987; 14(3): 221-233."),
      p("Borgan Ø, LiestØl K. A note on confidence intervals and bands for the survival function based on transformations."),
      p("Scandinavian Journal of Statistics 1990; 17(1): 35-41."),
      p("Nagashima K, Noma H, Sato Y, Gosho M. Sample size calculations for single-arm survival studies using transformations of the Kaplan-Meier estimator. Pharmaceutical Statistics 2020. DOI: 10.1002/pst.2090. [arXiv:2012.03355].")
    )
  )
}
