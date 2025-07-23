library(shiny)
library(bslib)


lapply(list.files("source", full.names = TRUE), source)

ui <- page_sidebar(
  layout = "sidebar-top",
  sidebar = sidebar(
    selectInput("tool", "Choose a function:", choices = c(
      "Test 1 Mean" = "one_mean_size",
      "Compare 2 Means" = "two_mean_size",
      "Compare k Means" = "k_mean_size",
      "Test 1 Proportion" = "one_prop_size",
      "Compare 2 Proportion" = "two_prop_size",
      "Compare Paired Proportion" = "pair_prop_size",
      "Compare k Proportion" = "k_prop_size",
      "Test Time-To-Event Data" = "coxph_size",
      "Test Odds Ratio" = "or_size",
      "SCCS, Alt-2" = "sccs_size",
      "One Normal" = "one_norm_size",
      "One Binomial" = "one_bino_size"
    )),
    uiOutput("param_ui")
  ),
  uiOutput("main_ui")
)

server <- function(input, output, session) {

  output$param_ui <- renderUI({
    switch(input$tool,
           "one_mean_size" = param_ui_one_mean_size(),
           "two_mean_size" = param_ui_two_mean_size(),
           "k_mean_size" = param_ui_k_mean_size(),
           "one_prop_size" = param_ui_one_prop_size(),
           "two_prop_size" = param_ui_two_prop_size(),
           "pair_prop_size" = param_ui_pair_prop_size(),
           "k_prop_size" = param_ui_k_prop_size(),
           "coxph_size" = param_ui_coxph_size(),
           "or_size" = param_ui_or_size(),
           "sccs_size" = param_ui_sccs_size(),
           "one_norm_size" = param_ui_one_norm_size(),
           "one_bino_size" = param_ui_one_bino_size()
    )
  })

  output$main_ui <- renderUI({
    switch(input$tool,
           "one_mean_size" = main_ui_one_mean_size(input$test_type),
           "two_mean_size" = main_ui_two_mean_size(input$test_type),
           "k_mean_size" = main_ui_k_mean_size(input$test_type),
           "one_prop_size" = main_ui_one_prop_size(input$test_type),
           "two_prop_size" = main_ui_two_prop_size(input$test_type),
           "pair_prop_size" = main_ui_pair_prop_size(input$test_type),
           "k_prop_size" = main_ui_k_prop_size(),
           "coxph_size" = main_ui_coxph_size(input$test_type),
           "or_size" = main_ui_or_size(input$test_type),
           "sccs_size" = main_ui_sccs_size(),
           "one_norm_size" = main_ui_one_norm_size(),
           "one_bino_size" = main_ui_one_bino_size()

    )
  })

  output$result <- renderText({
    res <- switch(input$tool,
                  "one_mean_size" = result_one_mean_size(input),
                  "two_mean_size" = result_two_mean_size(input),
                  "k_mean_size" = result_k_mean_size(input),
                  "one_prop_size" = result_one_prop_size(input),
                  "two_prop_size" = result_two_prop_size(input),
                  "pair_prop_size" = result_pair_prop_size(input),
                  "k_prop_size" = result_k_prop_size(input),
                  "coxph_size" = result_coxph_size(input),
                  "or_size" = result_or_size(input),
                  "sccs_size" = result_sccs_size(input),
                  "one_norm_size" = result_one_norm_size(input),
                  "one_bino_size" = result_one_bino_size(input)
    )
    if (input$mode == "size") {
      paste("Sample size: ", res)
    } else if (input$mode == "power") {
      paste("Power: ", res)
    }
  })
}

shinyApp(ui, server)
