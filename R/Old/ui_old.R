
library(shiny)
library(shinyFiles)

ui <- fluidPage(
  useShinyjs(),
  
  fluidRow(
    column(12,
           fluidRow(
             column(3,
                    id = "view_ongraph",
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
                                selected = 3)
             ),
             column(9,
                    fluidRow(
                      column(9,
                             sidebarPanel( width = 9,
                                           fluidRow(
                                             column(4, align = "center",
                                                    p(strong("Load Adjacency Matrix :"))),
                                             column(8, 
                                                    fileInput(inputId = "data_adjacency",NULL, accept = ".RData"))
                                           ),
                                           fluidRow(
                                             column(4, align = "center",
                                                    p(strong("Load Node Coordinates"))),
                                             column(8, 
                                                    fileInput(inputId = "data_coordinates", NULL, accept = ".RData"))
                                           ),
                                           fluidRow(
                                             column(4, align = "center",
                                                    p(strong("Load Event Coordinates"))),
                                             column(8, 
                                                    fileInput(inputId = "data_events", NULL, accept = ".RData"))
                                           )
                             ),
                             column(3,
                                    radioButtons(inputId = "select_net",
                                                 label= "Select network type",
                                                 choices = list("Undirected" = 1,
                                                                "Directed" = 2,
                                                                "Mixed" = 3),
                                                 selected = 1,
                                                 inline = TRUE),
                                    actionButton("load_net", "Load")
                             )
                      )
                    )
             )
           ),
           fluidRow(
             column(3,
                    h3("Node Info:"),
                    htmlOutput("node_info"),
                    h3("Edge Info:"),
                    htmlOutput("edge_info")
                    #dataTableOutput("nodes_data_from_shiny"),
                    #verbatimTextOutput("node_id"),
                    #verbatimTextOutput("node_int"),
             ),
             column(9,
                    visNetworkOutput("network")
             )
           )
    )
  )
)


