generateHeatmapUI <- function (id) {
  
  ns <- NS(id)
  
  tagList(
    h3("Heatmap"),
    fluidRow(    
      column(3,
             checkboxInput(ns("showSig"), "Show all significant features."),
             numericInput(ns("featureNumber"), "Select top n features for heatmap:", 50),
             selectInput(
               ns("heatScale"), 
               "Scale data by:", 
               choices = c("row", "column", "none"),
               selected = "none"
               ),
             ), 
      
      column(4, offset = 1,
             checkboxInput(ns("clustCol"), "Cluster Columns", TRUE),
             checkboxInput(ns("clustRow"), "Cluster Rows", TRUE),
      ),
      column(4,
             checkboxInput(ns("namesCol"), "Show Column Names", TRUE),
             checkboxInput(ns("namesRow"), "Show Row Names", TRUE),
      ),

    ),
    plotOutput(outputId = ns("heatmap")),
  )
  
}


generateHeatmapServer <- function(id, groupIdentities, transdf, ptable) {
  
  moduleServer(id, function(input, output, session) {
    
    plotHeight <- min(1000, 0.1 * length(rownames(workingDf)))
    ptable <- as.data.frame(ptable)
    
    output$heatmap <- renderPlot({
        
        ptable <- ptable[order(ptable[, 2]),]
        
        if (input$showSig == FALSE) {
          
          ptable <- ptable[1:input$featureNumber, 1] %>% sapply(as.character)
          
        } else {
          
          ptable <- ptable[which(ptable[, 2] < 0.05), 1] %>% sapply(as.character)
          
        }
        
        
        #create a df with the top 100 metabolites between mean Control and AD
        heatmap_data <- transdf[,c(which(as.character(names(transdf)) %in% ptable == TRUE | as.character(names(transdf)) == "Group"))]
        
        rowAnnot <- data.frame(Group = heatmap_data$Group)
        heatmap_data$Group <- NULL
        rownames(rowAnnot) <- rownames(heatmap_data)
        levels(rowAnnot$Group) <- levels(groupIdentities$Group)
        
        heatmap_data[heatmap_data == 0] <-1
        
        #Log transform the data
        heatmap_data <- heatmap_data %>% log2()
        
        #mean subtraction 
        heatmap_data <- heatmap_data - rowMeans(heatmap_data)
        
        tomerge <- data.frame(key = names(heatmap_data))
        tomerge <- merge(tomerge, key, by = "key")
        names(heatmap_data) <- tomerge$m.z
        
        
        heatmap_data %>%
          pheatmap(annotation_row = rowAnnot,
                   #annotation_colors = annotcolors,
                   scale = input$heatScale,
                   show_colnames = input$namesCol,
                   show_rownames = input$namesRow,
                   cluster_rows = input$clustRow,
                   cluster_cols = input$clustCol,
          )

    }, height = plotHeight)
    
  })
  
  
}