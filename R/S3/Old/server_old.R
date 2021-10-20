library(shiny)
library(shinyFiles) # Provides functionality for client-side navigation of the server side file system in shiny apps. 
library(shinyjs) # Perform common useful JavaScript operations in Shiny 
library(shiny.router) # It is a simple router for your Shiny apps
source("./main.R", local=TRUE)

server <- function(input, output, session) {
  output$show_file_panel <- reactive(TRUE)
  output$show_load_panel <- reactive(FALSE)
  output$show_net_panel <- reactive(FALSE)
  
  outputOptions(output, "show_file_panel", suspendWhenHidden = FALSE)
  outputOptions(output, "show_load_panel", suspendWhenHidden = FALSE)
  outputOptions(output, "show_net_panel", suspendWhenHidden = FALSE)
  
  
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
      
      # Adjacency matrix (undirected): Segmenting locations of the traffic network treated as the vertex set of the network.
      load(adj_mtx_path)
      print(adj_mtx_path)
      
      # Node coordinates: Georeferenced coordinates from 'castellon' nodes
      load(node_coords_path)
      
      # Event (crime coordinates)
      load(events_path)
    } 
    list(adj_mtx_path = adj_mtx_path, node_coords_path = node_coords_path, events_path=events_path)
  })
  
  
  observeEvent(input$load_net, {
    g <- NULL
    file_paths <- load_files()
    if(is.null(file_paths$adj_mtx_path) || is.null(file_paths$adj_mtx_path) || is.null(file_paths$events_path)){
      showModal(modalDialog(
        title = "Warning",
        "Please, select the files to create the network.",
        easyClose = TRUE,
        footer = NULL
      ))
    }else{
      output$show_file_panel <- reactive(FALSE)
      output$show_load_panel <- reactive(TRUE)
      output$show_net_panel <- reactive(FALSE)
      
      output$network <- renderVisNetwork({
        
        #subset of events
        crim <- crimes[11:111,] # From crimes, take 11 to 111 (both included)
        
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
        cat(paste(vertex_attr_names(intnet$graph), "\n"))
        #cat(paste(vertex_attr(intnet$graph)$xcoord[1], "\n"))
        g <<- intnet$graph
        edge_ids <- get.edge.ids(g, as.vector(t(get.edgelist(g)))) 
        
        nodes <- data.frame(id = paste(vertex_attr(g)$name),
                            x = vertex_attr(g)$xcoord,
                            y = vertex_attr(g)$ycoord,
                            label = paste(vertex_attr(g)$name))
        
        edges <- data.frame(id = edge_ids,
                            from = get.edgelist(g)[,1],
                            to = get.edgelist(g)[,2])
        
        
        visNetwork(nodes, edges)  %>%
          visIgraphLayout() %>%
          visEvents(selectNode = "function(nodes) {
                      Shiny.setInputValue('current_node_id', nodes.nodes);
                      ;}",
                    select = "function(edges) {
                      Shiny.setInputValue('current_edge_id', edges.edges);
                      ;}")
      })
      
      # observeEvent(input$current_node_id, {
      #   visNetworkProxy("network") %>%
      #     visGetNodes() 
      # })
      # 
      # observeEvent(input$current_edge_id, {
      #   visNetworkProxy("network") %>%
      #     visGetEdges()
      # })
      # 
      # node_model <- eventReactive(input$current_node_id, {
      #   if (!is.null(input$current_node_id)){
      #     node_int = vertex_attr(g, 'intensity', input$current_node_id)
      #     moran_i = vertex_attr(g, 'moran_i', input$current_node_id)
      #     getis_g = vertex_attr(g, 'getis_g', input$current_node_id)
      #     
      #     # return all object as a list
      #     list(node_id = input$current_node_id, node_intensity = node_int, node_moran = moran_i, node_g = getis_g)
      #   }
      # })
      # 
      # edge_model <- eventReactive(input$current_edge_id, {
      #   if (!is.null(input$current_edge_id)){
      #     edge_int = edge_attr(g, 'intensity', input$current_edge_id)
      #     
      #     # return all object as a list
      #     list(edge_id = input$current_edge_id, edge_intensity = edge_int)
      #   }
      # })
      # 
      # output$node_info <- renderText({
      #   paste0("<B>Id: </B>", node_model()$node_id, "<br>",
      #          "<B>Intensity: </B>",  node_model()$node_intensity, "<br>",
      #          "<B>Moran-i: </B>",  node_model()$node_moran, "<br>",
      #          "<B>Getis-g: </B>",  node_model()$node_g, "<br>")
      # }) 
      # 
      # output$edge_info <- renderText({
      #   paste0("<B>Id: </B>", input$current_edge_id, "<br>",
      #          "<B>Intensity: </B>",  vertex_attr(g, 'intensity', input$current_edge_id), "<br><br>")
      # })  
      # 
      # observe({
      #   if (!is.null(input$node_display)){
      #     if(input$node_display == 1){
      #       nodes <- data.frame(id = paste(vertex_attr(g)$name),
      #                           label = paste(vertex_attr(g)$name))
      #     }else if(input$node_display == 2){
      #       nodes <- data.frame(id = paste(vertex_attr(g)$name),
      #                           label = paste(round(vertex_attr(g)$intensity, 4)))
      #     }else if(input$node_display == 3){
      #       nodes <- data.frame(id = paste(vertex_attr(g)$name),
      #                           label = paste(round(vertex_attr(g)$moran_i, 4)))
      #     }else if(input$node_display == 4){
      #       nodes <- data.frame(id = paste(vertex_attr(g)$name),
      #                           label = paste(round(vertex_attr(g)$getis_g, 4)))
      #     }else{
      #       nodes <- data.frame(id = paste(vertex_attr(g)$name),
      #                           label = "")
      #     }
      #     visNetworkProxy("network") %>% visUpdateNodes(nodes)    
      #   }
      # })
      # 
      # observe({
      #   if (!is.null(input$edge_display)){
      #     edge_ids <- get.edge.ids(g, as.vector(t(get.edgelist(g)))) 
      #     
      #     if(input$edge_display == 1){
      #       edges <- data.frame(id = edge_ids,
      #                           from = get.edgelist(g)[,1],
      #                           to = get.edgelist(g)[,2],
      #                           label = paste(edge_ids))
      #     }else if(input$edge_display == 2){
      #       edges <- data.frame(id = edge_ids,
      #                           from = get.edgelist(g)[,1],
      #                           to = get.edgelist(g)[,2],
      #                           label = paste(round(edge_attr(g)$intensity, 4)))
      #     }else{
      #       edges <- data.frame(id = edge_ids,
      #                           from = get.edgelist(g)[,1],
      #                           to = get.edgelist(g)[,2],
      #                           label = " ")
      #     }
      #     visNetworkProxy("network") %>% visUpdateEdges(edges)    
      #   }
      # })
      # 
      # observe({
      #   print(input$network_selectedNodes)
      # })
    }
  })#End observeEvent load_net
  
  
  
  
  
}