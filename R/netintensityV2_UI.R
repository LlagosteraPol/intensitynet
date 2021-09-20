library(shiny)
library(shinyFiles)
library(shinythemes)
source("netintensityV2_UndirectedFactory.R")

# https://rstudio.github.io/shinythemes/ # Catalog of shiny themes

# Define UI
ui <- fluidPage(theme = shinytheme("cerulean"),
                navbarPage(
                  # theme = "cerulean",  # <--- To use a theme, uncomment this
                  "Intensitynet",
                  tabPanel("Navbar 1",
                           sidebarPanel(
                             tags$h3("Input:"),
                             shinyDirButton("dir", "Input Adjacency Matrix", "Upload"),
                             verbatimTextOutput("mtx_path", placeholder = TRUE),
                             textInput("txt1", "Given Name:", ""),
                             textInput("txt2", "Surname:", ""),
                             
                           ), # sidebarPanel
                           mainPanel(
                             h1("Header 1"),
                             
                             h4("Output 1"),
                             verbatimTextOutput("txtout"),
                             
                           ) # mainPanel
                           
                  ), # Navbar 1, tabPanel
                  tabPanel("Navbar 2", "This panel is intentionally left blank"),
                  tabPanel("Navbar 3", "This panel is intentionally left blank")
                  
                ) # navbarPage
) # fluidPage


# Define server function  
server <- function(input, output) {
  
  # load("Data/Castellon.RData")
  # 
  # undirected <- UndirectedFactory()
  # undirected_general <- undirected$constructGraph(Castellon, "general", nodes$cx, nodes$cy, crim$X, crim$Y)
  shinyDirChoose(
    input,
    'dir',
    roots = c(home = dirname(getwd())),
    filetypes = c('', 'txt', 'bigWig', "tsv", "csv", "bw", "RData")
  )
  
  global <- reactiveValues(datapath = dirname(getwd()))
  
  dir <- reactive(input$dir)
  
  output$mtx_path <- renderText({
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
  
  output$txtout <- renderText({
    paste( input$txt1, input$txt2, sep = " " )
  })
} # server


# Create Shiny object
shinyApp(ui = ui, server = server)