library(shiny)
library(bslib)


lapply(list.files("source", full.names = TRUE), source)

ui <- page_sidebar(
  sidebar = sidebar(
    width = "400px",
    open = list(
      desktop = "open",
      mobile = "always-above"
    ),
    selectInput("tool", "Test type:", choices = c(
      "Non-Inferiority" = "Non-Inferiority",
      "Superiority" = "Superiority",
      "One-sample" = "One-sample"
    )),
    uiOutput("param_ui")
  ),
  uiOutput("main_ui")
)

server <- function(input, output, session) {

  output$param_ui <- renderUI({
    switch(input$tool,
           "Non-Inferiority" = param_ui_twoSurvSampleSizeNI(),
           "Superiority" = param_ui_lakatosSampleSize(),
           "One-sample" = param_ui_oneSurvSampleSize()
    )
  })

  output$main_ui <- renderUI({
    switch(input$tool,
           "Non-Inferiority" = main_ui_twoSurvSampleSizeNI(),
           "Superiority" = main_ui_lakatosSampleSize(),
           "One-sample" = main_ui_oneSurvSampleSize()
    )
  })

  output$result <- renderTable({
    res <- switch(input$tool,
                  "Non-Inferiority" = result_twoSurvSampleSizeNI(input),
                  "Superiority" = result_lakatosSampleSize(input),
                  "One-sample" = result_oneSurvSampleSize(input)
    )
    Metric <- gsub("_", " ", names(res))
    data.frame(
      Metric,
      Value = unname(unlist(res)),
      row.names = NULL
    )
  }, colnames = FALSE)
}

shinyApp(ui, server)
