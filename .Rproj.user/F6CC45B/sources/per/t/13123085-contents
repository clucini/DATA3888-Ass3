shinyUI(fluidPage(
  titlePanel("Accuracy of KNN"),
  sidebarLayout(
    sidebarPanel(
      sliderInput("slider1", label = h3("Slider for value of K"), min = 1, 
                  max = 15, value = 5),
      actionButton(inputId = "button",
                   label = "Run experiment"),
      actionButton(inputId = "cum_button",
                   label = "Reset cumulative plot")
    ),
    mainPanel(
      shiny::plotOutput(outputId = "wind_plot"),
      shiny::plotOutput(outputId = "cum_plot")
    )
  )
))
