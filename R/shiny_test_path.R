library(shiny)
library(shinyFiles)

ui <- fluidPage( # Application title
  mainPanel(
    shinyDirButton("dir", "Input directory", "Upload"),
    verbatimTextOutput("dir", placeholder = TRUE)  
  ))

server <- function(input, output) {
  shinyDirChoose(
    input,
    'dir',
    roots = c(home = dirname(getwd())),
    filetypes = c('', 'txt', 'bigWig', "tsv", "csv", "bw", "RData")
  )
  
  global <- reactiveValues(datapath = dirname(getwd()))
  
  dir <- reactive(input$dir)
  
  output$dir <- renderText({
    normalizePath(global$datapath)
  })
  
  observeEvent(ignoreNULL = TRUE,
               eventExpr = {
                 input$dir
               },
               handlerExpr = {
                 if (!"path" %in% names(dir())){
                   return()
                 } 
                 home <- dirname(getwd())
                 global$datapath <-
                   file.path(home, paste(unlist(dir()$path[-1]), collapse = .Platform$file.sep))
               })
}

# Run the application
shinyApp(ui = ui, server = server)