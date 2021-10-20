
library(shiny)
library(shinyFiles)

ui <- fluidPage(
  uiOutput("file_panel"),
  uiOutput("load_panel"),
  uiOutput("net_panel")
)