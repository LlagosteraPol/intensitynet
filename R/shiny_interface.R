

library(shiny)
library(shinyFiles)
source("./main.R", local=TRUE)



ui <- fluidPage(
  fluidRow(
    column(3,
           selectInput(inputId = "select_net", 
                       h3("Select network type"), 
                       choices = list("Undirected" = 1, 
                                      "Directed" = 2,
                                      "Mixed" = 3), 
                       selected = 1))
  ),
  fluidPage(
    visNetworkOutput("network")
  ),
  column(3, 
         radioButtons("node_display", 
                      h3("Node display"), 
                      choices = c("id" = 1, 
                                  "intensity" = 2, 
                                  "moran i" = 3,
                                  "getis g" = 4,
                                  "none" = 5),
                      inline = TRUE)),
  
)

server <- function(input, output) {
  library(shiny)
  library(shinyFiles)
  source("./main.R", local=TRUE)
  # Adjacency matrix (undirected): Segmenting locations of the traffic network treated as the vertex set of the network.
  load("../../Data/Castellon.RData")
  
  # Node coordinates: Georeferenced coordinates from 'castellon' nodes
  load("../../Data/nodes.RData")
  
  # Event (crime coordinates)
  load("../../Data/crimes.RData")
  
  #subset of events
  crim <- crimes[11:111,] # From crimes, take 11 to 111 (both included)

  
  output$network <- renderVisNetwork({
    if(input$select_net == 1){
      intnet <- intensitynet(Castellon, nodes, crim)
    }else{
      Castellon_obj <-list(mtx = Castellon)
      class(Castellon_obj) <- "netTools"
      dirCastellon <-  Undirected2RandomDirectedAdjMtx(Castellon_obj)
      
      if(input$select_net == 2){
        intnet <- intensitynet(dirCastellon, nodes, crim, graph_type='directed')
      }else{
        intnet <- intensitynet(dirCastellon, nodes, crim, graph_type='mixed')
      }
    }
    
    intnet <- CalculateEventIntensities(intnet)
    intnet <- NodeLocalCorrelation(intnet, 'moran')
    intnet <- NodeLocalCorrelation(intnet, 'g')
    g <- intnet$graph
    
    nodes <- data.frame(id = paste(vertex_attr(g)$name),
                        x = vertex_attr(g)$xcoord,
                        y = vertex_attr(g)$ycoord)
    
    edges <- data.frame(from = get.edgelist(g)[,1],
                        to = get.edgelist(g)[,2])
    
    visNetwork(nodes, edges)  %>%
      visIgraphLayout()
    })
  
  observe({
    if (!is.null(input$node_display)){
      if(input$node_display == 1){
        nodes <- data.frame(id = paste(vertex_attr(g)$name),
                            label = paste(vertex_attr(g)$name))
      }else if(input$node_display == 2){
        nodes <- data.frame(id = paste(vertex_attr(g)$name),
                            label = paste(round(vertex_attr(g)$intensity, 4)))
      }else if(input$node_display == 3){
        nodes <- data.frame(id = paste(vertex_attr(g)$name),
                            label = paste(round(vertex_attr(g)$moran_i, 4)))
      }else if(input$node_display == 4){
        nodes <- data.frame(id = paste(vertex_attr(g)$name),
                            label = paste(round(vertex_attr(g)$getis_g, 4)))
      }else{
        nodes <- data.frame(id = paste(vertex_attr(g)$name),
                            label = "")
      }
      visNetworkProxy("network") %>% visUpdateNodes(nodes)    
    }
  })
  
}

# Run the application
shinyApp(ui = ui, 
         server = server
         )