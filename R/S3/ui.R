
library(shiny)
library(shinyFiles)

ui <- fluidPage(
  useShinyjs(),
  
  conditionalPanel(condition = "output.show_file_panel",
                   h1("Choose files"),
                   sidebarPanel( width = 12,
                                 fileInput(inputId = "data_adjacency",label = "Load Adjacency Matrix :", accept = ".RData"),
                                 fileInput(inputId = "data_coordinates", label = "Load Node Coordinates", accept = ".RData"),
                                 fileInput(inputId = "data_events", label = "Load Event Coordinates", accept = ".RData")
                   ),
                   radioButtons(inputId = "select_net",
                                label= "Select network type",
                                choices = list("Undirected" = 1,
                                               "Directed" = 2,
                                               "Mixed" = 3),
                                selected = 1,
                                inline = TRUE),
                   actionButton("load_net", "Load")
  ),
  
  conditionalPanel(condition = "output.show_load_panel",
                   h1("LOADING..."),
                   actionButton("proceed", "proceed")
  ),
  
  conditionalPanel(condition = "output.show_net_panel",
                   fluidRow(
                     column(3,
                            selectInput(inputId = "node_display", 
                                        label = "Node display:",
                                        choices = c("id" = 1, 
                                                    "intensity" = 2, 
                                                    "moran i" = 3,
                                                    "getis g" = 4,
                                                    "none" = 5),
                                        selected = 1),
                            selectInput(inputId = "edge_display", 
                                        label = "Edge display", 
                                        choices = c("id" = 1, 
                                                    "intensity" = 2, 
                                                    "none" = 3),
                                        selected = 3),
                            h3("Node Info:"),
                            htmlOutput("node_info"),
                            h3("Edge Info:"),
                            htmlOutput("edge_info")
                            #dataTableOutput("nodes_data_from_shiny"),
                            #verbatimTextOutput("node_id"),
                            #verbatimTextOutput("node_int"),
                     ),
                     column(9,
                            visNetworkOutput("network",  height = "100vh")
                     )
                   )
  )
)


