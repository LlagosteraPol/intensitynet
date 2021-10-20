library(shiny)
library(shinyFiles) # Provides functionality for client-side navigation of the server side file system in shiny apps. 
library(shinyjs) # Perform common useful JavaScript operations in Shiny 
#library(shiny.router) # It is a simple router for your Shiny apps
source("./main.R", local=TRUE)


gen_net <- function(adj_mtx, node_coords, events){
  intnet <- intensitynet(adj_mtx, node_coords, events)
  intnet <- CalculateEventIntensities(intnet)
  intnet <- NodeLocalCorrelation(intnet, 'moran')
  intnet <- NodeLocalCorrelation(intnet, 'g')
  cat(paste(vertex_attr_names(intnet$graph), "\n"))
  #cat(paste(vertex_attr(intnet$graph)$xcoord[1], "\n"))
  g <- intnet$graph
  edge_ids <- get.edge.ids(g, as.vector(t(get.edgelist(g)))) 
  
  nodes <- data.frame(id = paste(vertex_attr(g)$name),
                      x = vertex_attr(g)$xcoord,
                      y = vertex_attr(g)$ycoord,
                      label = paste(vertex_attr(g)$name))
  
  edges <- data.frame(id = edge_ids,
                      from = get.edgelist(g)[,1],
                      to = get.edgelist(g)[,2])
  
  list(nodes = nodes, edges = edges, graph = g)
}

server <- function(input, output, session) {
  g <- NULL
  show_panel <- "file_panel"
  
  output$file_panel <- renderUI({
    if(show_panel == "file_panel"){
      div(
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
      )
    } else return(div())
  })
  
  output$load_panel <- renderUI({
    if(input$load_net == 0){
      return(div())
      
    }else {
      output$file_panel <- renderUI({div()})
      div(h1("LOADING..."))
    }
  })
  
  output$net_panel <- renderUI({
    if(show_panel == "net_panel"){
      output$load_panel <- renderUI({div()})
      div(
        fluidRow(
          column(12,
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
    }else return(div())
  })
  
  
  load_files <- reactive({
    adj_mtx_path <- input$data_adjacency
    node_coords_path <- input$data_coordinates
    events_path <- input$data_events
    
    if (is.null(adj_mtx_path) || is.null(adj_mtx_path) || is.null(events_path)){
      #Default values
      adj_mtx_path <- "../../Data/Castellon.RData"
      node_coords_path <- "../../Data/nodes.RData"
      events_path <- "../../Data/crimes.RData"
      #return(NULL) # Uncomment if not default values
    } 
    list(adj_mtx_path = adj_mtx_path, node_coords_path = node_coords_path, events_path=events_path)
  })
  
  observeEvent(input$load_net, {
    
    file_paths <- load_files()
    if(is.null(file_paths$adj_mtx_path) || is.null(file_paths$adj_mtx_path) || is.null(file_paths$events_path)){
      showModal(modalDialog(
        title = "Warning",
        "Please, select the files to create the network.",
        easyClose = TRUE,
        footer = NULL
      ))
    }else{
      show_panel <<- "net_panel"
      
      # Adjacency matrix (undirected): Segmenting locations of the traffic network treated as the vertex set of the network.
      load(file_paths$adj_mtx_path)
      print(file_paths$adj_mtx_path)
      
      # Node coordinates: Georeferenced coordinates from 'castellon' nodes
      cc <- load(file_paths$node_coords_path)
      
      # Event (crime coordinates)
      load(file_paths$events_path)
      #subset of events
      crim <- crimes[11:111,] # From crimes, take 11 to 111 (both included)
      
      net_data <- gen_net(Castellon, nodes, crim)
      nodes <- net_data$nodes
      edges <- net_data$edges
      g <- net_data$graph
      
      # if(input$select_net == 1){
      #   intnet <- intensitynet(Castellon, nodes, crim)
      # }else{
      #   Castellon_obj <-list(mtx = Castellon)
      #   class(Castellon_obj) <- "netTools"
      #   dirCastellon <-  Undirected2RandomDirectedAdjMtx(Castellon_obj)
      #   
      #   if(input$select_net == 2){
      #     intnet <- intensitynet(dirCastellon, nodes, crim, graph_type='directed')
      #   }else{
      #     intnet <- intensitynet(dirCastellon, nodes, crim, graph_type='mixed')
      #   }
      # }
      # output$show_file_panel <- reactive(FALSE)
      # output$show_load_panel <- reactive(FALSE)
      # output$show_net_panel <- reactive(TRUE)
      # outputOptions(output, "show_file_panel", suspendWhenHidden = FALSE)
      # outputOptions(output, "show_load_panel", suspendWhenHidden = FALSE)
      # outputOptions(output, "show_net_panel", suspendWhenHidden = FALSE)
      
    } # End else
  }) #End observeEvent load_net
  
  observeEvent(g,{
    output$show_file_panel <- reactive(FALSE)
    output$show_load_panel <- reactive(FALSE)
    output$show_net_panel <- reactive(TRUE)
    outputOptions(output, "show_file_panel", suspendWhenHidden = FALSE)
    outputOptions(output, "show_load_panel", suspendWhenHidden = FALSE)
    outputOptions(output, "show_net_panel", suspendWhenHidden = FALSE)
    
    output$network <- renderVisNetwork({
      visNetwork(nodes, edges)  %>%
        visIgraphLayout() %>%
        visEvents(selectNode = "function(nodes) {
                      Shiny.setInputValue('current_node_id', nodes.nodes);
                      ;}",
                  select = "function(edges) {
                      Shiny.setInputValue('current_edge_id', edges.edges);
                      ;}")
    })
    
    observeEvent(input$current_node_id, {
      visNetworkProxy("network") %>%
        visGetNodes()
    })
    
    observeEvent(input$current_edge_id, {
      visNetworkProxy("network") %>%
        visGetEdges()
    })
    
    node_model <- eventReactive(input$current_node_id, {
      if (!is.null(input$current_node_id)){
        node_int = vertex_attr(g, 'intensity', input$current_node_id)
        moran_i = vertex_attr(g, 'moran_i', input$current_node_id)
        getis_g = vertex_attr(g, 'getis_g', input$current_node_id)
        
        # return all object as a list
        list(node_id = input$current_node_id, node_intensity = node_int, node_moran = moran_i, node_g = getis_g)
      }
    })
    
    edge_model <- eventReactive(input$current_edge_id, {
      if (!is.null(input$current_edge_id)){
        edge_int = edge_attr(g, 'intensity', input$current_edge_id)
        
        # return all object as a list
        list(edge_id = input$current_edge_id, edge_intensity = edge_int)
      }
    })
    
    output$node_info <- renderText({
      paste0("<B>Id: </B>", node_model()$node_id, "<br>",
             "<B>Intensity: </B>",  node_model()$node_intensity, "<br>",
             "<B>Moran-i: </B>",  node_model()$node_moran, "<br>",
             "<B>Getis-g: </B>",  node_model()$node_g, "<br>")
    })
    
    output$edge_info <- renderText({
      paste0("<B>Id: </B>", input$current_edge_id, "<br>",
             "<B>Intensity: </B>",  vertex_attr(g, 'intensity', input$current_edge_id), "<br><br>")
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
    
    observe({
      if (!is.null(input$edge_display)){
        edge_ids <- get.edge.ids(g, as.vector(t(get.edgelist(g))))
        
        if(input$edge_display == 1){
          edges <- data.frame(id = edge_ids,
                              from = get.edgelist(g)[,1],
                              to = get.edgelist(g)[,2],
                              label = paste(edge_ids))
        }else if(input$edge_display == 2){
          edges <- data.frame(id = edge_ids,
                              from = get.edgelist(g)[,1],
                              to = get.edgelist(g)[,2],
                              label = paste(round(edge_attr(g)$intensity, 4)))
        }else{
          edges <- data.frame(id = edge_ids,
                              from = get.edgelist(g)[,1],
                              to = get.edgelist(g)[,2],
                              label = " ")
        }
        visNetworkProxy("network") %>% visUpdateEdges(edges)
      }
    })
    
    observe({
      print(input$network_selectedNodes)
    })
  })
}