generateLimmaUI <- function (id) {
  
  ns <- NS(id)
  
  tagList(
    h3("Limma & Volcano Plot"),
    tableOutput(ns("toptable")),
    downloadButton(ns('downloadToptable'), "Download Top Table"),
    plotOutput(ns("volcanoPlot")),
  )
  
}


generateLimmaServer <- function(id, transdf, plotTitle) {
  
  moduleServer(id, function(input, output, session) {
    
    
    top.table <- isolate({
        
        req(transdf)
        
        #Expression data
        assayData <- as.matrix(t(transdf[,1:(length(names(transdf))-1)]))
        
        #phenotype data
        pData <- data.frame(Group = transdf$Group)
        pData <- new("AnnotatedDataFrame", data = pData)
        sampleNames(pData) <- colnames(assayData)
        
        #Feature Data
        FeatureData <- as.data.frame(names(transdf[, 1:(length(transdf)-1)]))
        names(FeatureData) <- "key"
        FeatureData <- new("AnnotatedDataFrame", data=FeatureData)
        featureNames(FeatureData) <- rownames(assayData)
        
        
        eSet <- ExpressionSet(assayData = assayData,
                              phenoData = pData,
                              featureData = FeatureData)
        
        
        #Create a model matrix
        design <- model.matrix(~Group, data=pData(eSet))
        
        # Fit the model
        fit <- lmFit(eSet, design)
        
        # Calculate the t-stat 
        fit <- eBayes(fit)
        
        # Summarize the results
        #results <- decideTests(fit[,"Groupgroup2"])
        
        top.table <- topTable(fit, sort.by = "P", n = Inf)
        top.table <- merge(key, top.table, by = "key")
        top.table$key <- NULL
        top.table <- top.table[order(top.table$P),]

    })
    
    output$toptable <- renderTable({
        
        req(top.table)
        head(top.table)
      
    })
    
    
    output$downloadToptable <- downloadHandler(
      filename = function(){paste(plotTitle, "- top table.csv")}, 
      content = function(fname){
        
        write.csv(top.table, fname, row.names = FALSE)
        
      }
    )
    
    # Generate Volcano ----
    
    output$volcanoPlot <- renderPlot({
        
        req(top.table)
        
        top.table$difffex <- "NO"
        top.table$difffex[top.table$logFC > 0.6 & top.table$adj.P.Val < 0.05] <- "UP"
        # if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
        top.table$difffex[top.table$logFC < -0.6 & top.table$adj.P.Val < 0.05] <- "DOWN"
        
        top.table$label[top.table$difffex != "NO"] <- top.table$m.z[top.table$difffex != "NO"]
        
        ggplot(data = top.table, aes(x = logFC, y = -log10(adj.P.Val), color = difffex, label = label)) +
          geom_point() +
          geom_label_repel(max.overlaps = 20) +
          theme_bw() +
          scale_color_manual(values = c("blue", "black", "red")) +
          geom_vline(xintercept = c(-0.6, 0.6), col = "red") +
          geom_hline(yintercept = -log10(0.05), col = "red") + 
          labs(x = "Log2 of Fold Change", 
               y = "Significance (-log10P)",
               color = "Enrichment")
        
      
    })
    
  })
}
