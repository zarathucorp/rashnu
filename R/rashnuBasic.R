#' Interactive Sample Size Calculator for Survival Studies (Shiny App)
#'
#' Launches a Shiny web application that calculates required sample sizes and expected event numbers for different types
#' of survival analysis designs:
#' \itemize{
#'   \item Two-group Non-Inferiority
#'   \item Two-group Superiority (Lakatos method)
#'   \item One-sample survival test (with transformation methods)
#' }
#'
#' Users can specify survival probabilities, accrual and follow-up durations, allocation ratios,
#' non-inferiority margins, transformation methods, and test types. The app dynamically adjusts input UI
#' based on the selected design and displays results in a data table format.
#'
#' @details
#' \strong{Test Types:}
#' \itemize{
#'   \item \code{"ni"} - Non-Inferiority (two-group exponential survival comparison)
#'   \item \code{"sup"} - Superiority (Lakatos method with logrank/Gehan/Tarone-Ware weighting)
#'   \item \code{"one"} - One-sample survival test with multiple transformation options
#' }
#'
#' \strong{Included References:}
#' \itemize{
#'   \item Jung SH, Chow SC. Journal of Biopharmaceutical Statistics, 2012.
#'   \item Lakatos E. Biometrics, 1988.
#'   \item Lakatos & Lan. Statistics in Medicine, 1992.
#'   \item Fleming & Harrington. Counting Processes and Survival Analysis, 1991.
#'   \item Borgan Ø, Andersen PK et al. Springer-Verlag, 1993.
#'   \item Nagashima et al. Pharmaceutical Statistics, 2020.
#' }
#'
#' @note Requires associated functions \code{twoSurvSampleSizeNI()}, \code{lakatosSampleSize()}, and \code{oneSurvSampleSize()}
#' to be defined in the environment. Assumes a CSS file is available at \code{"www/style.css"} for custom styling.
#'
#' @return Launches a Shiny app in the default browser.
#'
#' @examples
#' if (interactive()) {
#'   rashnuBasic()
#' }
#'
#' @import shiny
#' @importFrom DT datatable DTOutput renderDT
#' @export
rashnuBasic <- function(){

  ui <- fluidPage(
    shiny::includeCSS(system.file("www/style.css", package = "rashnu")),
    titlePanel("Sample Size Calculator"),
    sidebarLayout(
      sidebarPanel(
        radioButtons("test_type", "Test Type:",
                     choices = c("Non-Inferiority" = "ni", "Superiority" = "sup", "One-sample" = "one"),
                     selected = "ni", inline = T),
        numericInput("syear", "Survival Time :", value = 12),
        uiOutput("surv_ui"),
        uiOutput("alloc_ui"),
        fluidRow(
          column(6,
                 numericInput("accrual", "Accrual Time :", value = 24)
          ),
          column(6,
                 numericInput("follow", "Follow-up Time :", value = 24)
          )
        ),
        uiOutput("alpha_ui"),
        numericInput("power", "Power (1 - Beta):", value = 0.8),
        conditionalPanel(
          condition = "input.test_type == 'ni'",
          numericInput("margin", "Non-inferiority Margin :", value = 1.3)
        ),
        conditionalPanel(
          condition = "input.test_type == 'sup'",
          selectInput("method", "Test Method:", choices = c("logrank", "gehan", "tarone-ware")),
          selectInput("side", "Hypothesis:", choices = c("two.sided", "one.sided"))
        ),
        conditionalPanel(
          condition = "input.test_type=='ni'",
          selectInput("side", "Hypothesis:", choices = "one.sided")
        ),
        conditionalPanel(
          condition = "input.test_type=='one'",
          selectInput("side", "Hypothesis:", choices = c("two.sided", "one.sided")),
          selectInput("trans", "Transformation:", choices = c("arcsin", "log-log", "logit", "log", "log-swog", "identity"))
        ),
        actionButton("calc", "Calculate")
      ),
      mainPanel(

         DTOutput("result_table"),
         tags$h3("Reference"),
         tags$div("Jung SH, Chow SC. On sample size calculation for comparing survival curves under general hypothesis testing. Journal of Biopharmaceutical Statistics 2012; 22(3):485–495."),

         tags$div(" Lakatos E. Sample sizes based on the log-rank statistic in complex clinical trials. Biometrics 1988; 44:229–241."),

         tags$div(" Lakatos E, Lan KK. A comparison of sample size methods for the logrank statistic. Statistics in Medicine 1992; 11(2):179–191."),

         tags$div(" Fleming TR, Harrington DP. Counting Processes and Survival Analysis. New York: Wiley, 1991, 236–237, Example 6.3.1."),

         tags$div(" Andersen PK, Borgan Ø, Gill RD, Keiding N. Statistical Models Based on Counting Processes. New York: Springer-Verlag, 1993, 176–287, Section IV.1–3."),

         tags$div(" Bie O, Borgan Ø, Liestøl K. Confidence intervals and confidence bands for the cumulative hazard rate function and their small sample properties. Scandinavian Journal of Statistics 1987; 14(3): 221–233."),

         tags$div(" Borgan Ø, Liestøl K. A note on confidence intervals and bands for the survival function based on transformations. Scandinavian Journal of Statistics 1990; 17(1): 35–41."),

         tags$div(" Nagashima K, Noma H, Sato Y, Gosho M. Sample size calculations for single-arm survival studies using transformations of the Kaplan–Meier estimator. Pharmaceutical Statistics 2020. DOI: 10.1002/pst.2090. [arXiv:2012.03355].")

      )
    )
  )

  server <- function(input, output) {
    observeEvent(input$calc, {


      res_df <- NULL

      if (input$test_type == "ni") {
        res <- tryCatch({
          twoSurvSampleSizeNI(
            syear = input$syear,
            yrsurv1 = input$yrsurv1,
            yrsurv2 = input$yrsurv2,
            accrualTime = input$accrual,
            followTime = input$follow,
            alloc = input$alloc,
            alpha = input$alpha,
            power = input$power,
            margin = input$margin
          )
        }, error = function(e) NULL)

        if (!is.null(res)) {

          res_df <- data.frame(
            Metric = names(res),
            Value = unname(unlist(res))
          )
        } else {
          res_df <- data.frame(Metric = "Error", Value = "Invalid or missing output from NI function")
        }

      } else if(input$test_type == "sup") {
        res <- lakatosSampleSize(
          syear = input$syear,
          yrsurv1 = input$yrsurv1,
          yrsurv2 = input$yrsurv2,
          accrualTime = input$accrual,
          followTime = input$follow,
          alloc = input$alloc,
          alpha = input$alpha,
          power = input$power,
          method = input$method,
          side = input$side
        )
        if (is.null(res$error)) {
          res_df <- data.frame(
            Metric = names(res),
            Value = unname(unlist(res))
          )
        } else {
          res_df <- data.frame(Metric = "Error", Value = res$error)
        }
      }else if(input$test_type =="one"){
        res <- tryCatch({oneSurvSampleSize(
          survTime = input$syear,
          p1 = input$p1,
          p2 = input$p2,
          accrualTime = input$accrual,
          followTime = input$follow,
          alpha = input$alpha,
          power = input$power,
          side = input$side,
          method = input$trans
        )}, error = function(e) NULL)

        if (!is.null(res)) {
          res_df <- data.frame(
            Metric = names(res),
            Value = unname(unlist(res))
          )
        } else {
          res_df <- data.frame(Metric = "Error", Value = "Invalid or missing output from NI function")
        }
      }

      output$result_table <- renderDT({
        datatable(
          res_df,
          colnames = c("Metric", "N"),
          options = list(
            dom = 't',
            ordering = FALSE,
            columnDefs = list(list(className = 'dt-center', targets = '_all'))
          ),
          rownames = FALSE
        )
      })
    })



    output$alpha_ui <- renderUI({
      if (input$test_type == "ni") {
        numericInput("alpha", "Significance Level (alpha):", value = 0.025)
      } else {
        numericInput("alpha", "Significance Level (alpha):", value = 0.05)
      }
    })

    output$alloc_ui <- renderUI({
      if (input$test_type %in% c("ni","sup")){
        numericInput("alloc", "Allocation Ratio :", value = 1)
      }
    })

    output$surv_ui <- renderUI({
      if (input$test_type %in% c("ni", "sup")){
        fluidRow(
          column(6,
                 numericInput("yrsurv1", "Survival Probability (Standard Group):", value = 0.305)
          ),
          column(6,
                 numericInput("yrsurv2", "Survival Probability (Test Group):", value = 0.435)
          )
        )
      }else{
        fluidRow(
          column(6,
                 numericInput("p1", "Null survival probability:", value = 0.305)
          ),
          column(6,
                 numericInput("p2", "Alternative survival probability", value = 0.435)
          )
        )
      }
    })

  }

  shinyApp(ui, server)

}
