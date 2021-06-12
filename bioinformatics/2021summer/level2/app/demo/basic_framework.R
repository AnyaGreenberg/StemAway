library(shiny)
library(ggplot2)

datasets <- c('economics', 'faithfuld', 'seals')

ui <- fluidPage(
  h1('Demo'),
  selectInput('dataset', label='Dataset', choices=datasets),
  verbatimTextOutput('summary'),
  tableOutput('table'),
  
  
  # Exercises 1.8
  h1('Exercises 1.8'),
  h2('#1'),
  textInput('name', "What's your name?"),
  textOutput('greeting'),
  
  h2('#2 and 3'),
  sliderInput('x', label = 'If x is', min=1, max=50, value=30),
  sliderInput('y', label='and y is', min=1, max=50, value=5),
  'then, x times y is',
  textOutput('product'),
  
  h2('#4'),
  'and (x*y)+5 is', textOutput('product_plus5'),
  'and (x*y)+10 is', textOutput('product_plus10'),
  
  h2('#5'),
  plotOutput('plot')
  
)

server <- function(input, output, session) {
  # Create reactive expression (function like functions and variables in Rshiny)
  dataset <- reactive({
    get(input$dataset, 'package:ggplot2')
  })
  
  # Paired with verbatimTextOutput()
  output$summary <- renderPrint({
    summary(dataset())
  })
  
  # Paired with tableOutput()
  output$table <- renderTable({
    head(dataset())
  })
  
  
  # Exercises 1.8
  output$greeting <- renderText({
    paste0('Hello ', input$name, '!')
  })
  
  times <- reactive({
    input$x*input$y
  })
  output$product <- renderText({
    times()
  })
  
  output$product_plus5 <- renderText({
    times()+5
  })
  
  output$product_plus10 <- renderText({
    times()+10
  })
  
  output$plot <- renderPlot({
    plot(dataset())
  }, res=96)
}

shinyApp(ui, server)
