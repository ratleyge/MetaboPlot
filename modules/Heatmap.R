generateHeatmapUI <- function (id) {
  
  ns <- NS(id)
  
  tagList(
    h3("Heatmap"),
    numericInput(ns("featureNumber"), "Select top n features for heatmap:", 50),
    plotOutput(outputId = ns("heatmap")),
  )
  
}


generateHeatmapServer <- function(id, groupIdentities, transdf, ptable) {
  
  moduleServer(id, function(input, output, session) {
    
    output$heatmap <- renderPlot({
        
        ptable <- ptable[order(ptable$pvalue),]
        #ptable <- ptable[which(ptable$pvalue < 0.05), 1] %>% sapply(as.character)
        ptable <- ptable[1:input$featureNumber, 1] %>% sapply(as.character)
        
        #create a df with the top 100 metabolites between mean Control and AD
        heatmap_data <- transdf[,c(which(as.character(names(transdf)) %in% ptable == TRUE | as.character(names(transdf)) == "Group"))]
        
        rowAnnot <- data.frame(Group = heatmap_data$Group)
        heatmap_data$Group <- NULL
        rownames(rowAnnot) <- rownames(heatmap_data)
        levels(rowAnnot$Group) <- levels(groupIdentities$Group)
        
        Group <- c("#00425A", "#FC7300")
        names(Group) <- c(levels(groupIdentities$Group)[1], levels(groupIdentities$Group)[2])
        annotcolors <- list(Group = Group)
        
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
                   annotation_colors = annotcolors,
                   scale = "column",
                   show_colnames = T,
                   show_rownames = F,
                   #cellheight = 20,
                   #height = 10,
                   #width = 10
          )
      
      
    })
    
  })
  
  
}