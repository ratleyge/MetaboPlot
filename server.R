server <- function(input, output, session) {
  
  options(shiny.maxRequestSize=100*1024^2)
  
  myData <- reactiveValues()
  
  # Outputs 
  output$viewInputTable <- renderTable({
    
    req(input$file1)
    workingDf <- read.csv(input$file1$datapath)
    
    workingDf$Bucket.label <- NULL
    workingDf$RT <- NULL
    workingDf$Name <- NULL
    workingDf$Formula <- NULL
    
    groupIdentities <- data.frame(workingDf[1, ])
    groupIdentities$m.z <- NULL
    groupIdentities <- as.data.frame(t(groupIdentities))
    names(groupIdentities) <- "Group"
    groupIdentities$Group <- as.factor(groupIdentities$Group)
    myData$groupIdentities <- groupIdentities
    
    workingDf <- workingDf[-1,]
    workingDf <<- workingDf
    
    head(workingDf[, 1:6])
    
  }) 
  
  
  output$Groups <- renderText({ 
    
    groupIdentities <- req(myData$groupIdentities)
    
    if (length(levels(groupIdentities$Group)) > 2) {
      
      paste0("<p>The data you uploaded has more than two groups: ", 
             paste(levels(groupIdentities$Group), collapse = ", "), 
             ". Please choose a different dataset.</p>")
      
    } else {
      
      paste("<p>Control:", levels(groupIdentities$Group)[1], "</p>", 
            "<p>Experimental:", levels(groupIdentities$Group)[2], "</p>")
      
    }
    
  })
  
  
  # Need to fix this so that people can relevel after nmds is generated and outliers are removed
  observeEvent(input$relevel, {
    
    if (is.null(myData$groupIdentities_outliersRemoved)) {
      
      groupIdentities <- myData$groupIdentities
      groupIdentities$Group <- relevel(groupIdentities$Group, levels(groupIdentities$Group)[2])
      myData$groupIdentities <- groupIdentities
      
    } else {  
      
      groupIdentities <- myData$groupIdentities_outliersRemoved
      groupIdentities$Group <- relevel(groupIdentities$Group, levels(groupIdentities$Group)[2])
      myData$groupIdentities_outliersRemoved <- groupIdentities
      
    }
    
  })
  
  
  observeEvent(input$NMDS, {
    
    output$NMDS <- renderPlot({
      
      groupIdentities <- req(myData$groupIdentities)
      
      showModal(modalDialog("Generating NMDS plot. This may take a moment.", footer=NULL))
      
      # Transpose Data and set column names as m/z
      transDf <- as.data.frame(t(workingDf))
      
      # Turn all the intensities into numerics and format as numeric matrix
      transDf <- sapply(transDf, as.numeric)
      colnames(transDf) <- workingDf$m.z
      transDf <- transDf[-1,]
      rownames(transDf) <- rownames(groupIdentities)
      
      # Run NMDS
      set.seed(123)
      nmds <- metaMDS(transDf, distance = "bray")
      nmds
      
      # Extract NMDS scores (x and y coordinates)
      data.scores <- as.data.frame(scores(nmds)$sites)
      
      if(input$outlierRemove == TRUE){
        
        # Identify outlier samples based on nmds values and safe a copy of them in directory
        outliers <- rbind(data.scores[which((data.scores$NMDS1 %in% boxplot.stats(data.scores$NMDS1)$out)), ],
                          data.scores[which((data.scores$NMDS2 %in% boxplot.stats(data.scores$NMDS2)$out)), ])
        # write.csv(outliers, paste(i, "- Outliers"), row.names = TRUE)
        
        # remove outliers from all data frames
        data.scores <- data.scores[which(!(rownames(data.scores) %in% rownames(outliers))), ]
        groupIdentities <- data.frame(Group = groupIdentities[which(!(rownames(groupIdentities) %in% rownames(outliers))), ])
        rownames(groupIdentities) <- rownames(data.scores)
        transDf <- transDf[which(!(rownames(transDf) %in% rownames(outliers))), ]
        workingDf <- workingDf[, which(!(colnames(workingDf) %in% rownames(outliers)))]
        
        myData$outliers <- rownames(outliers)
        workingDf <<- workingDf
        myData$groupIdentities_outliersRemoved  <- groupIdentities

      }
      
      
      # Create NMDS plot ###########################
      
      # Add group identities 
      data.scores$Group <- groupIdentities[, 1]
      head(data.scores)
      
      # Calculate ANOSIM statistic (will take a few min)
      
      
      # Now we can plot our NMDS in ggplot2
      xx <- ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) + 
        geom_point(size = 3, aes( shape = Group, colour = Group)) + 
        stat_ellipse(aes(x=NMDS1, y=NMDS2,color=Group),level = 0.95) +
        theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
              axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
              legend.text = element_text(size = 12, face ="bold", colour ="black"), 
              legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
              axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
              legend.title = element_text(size = 14, colour = "black", face = "bold"), 
              panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2),
              plot.title = element_text(face = "bold", size = 20, colour = "black"),
              legend.key=element_blank()) + 
        labs(x = "NMDS1", colour = "Group", y = "NMDS2", shape = "Group", title = paste(input$plotTitles, "NMDS")) 
      
      
      if (input$ANOSIM == TRUE) {
        
        ano <- anosim(transDf, groupIdentities[, 1], distance = "bray", permutations = 9999)
        xx <- xx + labs(subtitle = paste("ANOSIM stat:", 
                                         signif(ano$statistic, digits = 2),
                                         "Significance:", 
                                         signif(ano$signif, digits = 2)))
        
      }
      
      
      removeModal()
      print(xx)
      
    }, height = 400)
    
    
    output$outliers <- renderText(
      
      if(input$outlierRemove == TRUE){
      paste(
        "<p>The following samples were idenetified as outliers and are not pictured in the NMDS plot:</p> <ul> <li>",
        paste(myData$outliers, collapse='</li><li>'),
        "</li></ul><p>These samples will be removed from further analysis.</p>"
      )
        }
    )
    
  })
  
  
  observeEvent(input$heatmapButton, {
    
    ptable <- reactive({
      
      req(workingDf)
      
      if (is.null(myData$groupIdentities_outliersRemoved)) {
        
        groupIdentities <- myData$groupIdentities
        
      } else {   groupIdentities <- myData$groupIdentities_outliersRemoved   }
      
      
      workingDf <- cbind(key = c(paste0("met",seq(1:(length(rownames(workingDf)))))), workingDf)
      key <<- workingDf[, 1:2]
      workingDf$m.z <- NULL
      
      transdf <- as.data.frame(t(workingDf[, c(2:length(colnames(workingDf)))]))
      transdf <- as.data.frame(sapply(transdf, as.numeric))
      colnames(transdf) <- workingDf$key
      rownames(transdf) <- rownames(groupIdentities)
      
      
      # Create an alias for Groups so that we can call it in the p-table generation
      groupIdentities$AutoGroup <- 0
      groupIdentities$Group <- as.factor(groupIdentities$Group)
      groupIdentities[which(groupIdentities$Group == levels(groupIdentities$Group)[1]),]$AutoGroup <- "groupOne"
      groupIdentities[which(groupIdentities$Group == levels(groupIdentities$Group)[2]),]$AutoGroup <- "groupTwo"
      transdf$Group <- groupIdentities$AutoGroup
      transdf$Group <- as.factor(transdf$Group)
      transdf <<- transdf
      
      #Create a table of pvalues 
      ptable <- transdf %>% #select only a few columns to start with
        gather(key = m.z, value = intensity, -Group) %>% #melt the data into a long format
        group_by(Group, m.z) %>% #group intensities by metabolite and cohort
        dplyr::summarize(intensity = list(intensity)) %>% #put the intensities for all samples in a cohort for a particular metabolite in one column
        spread(Group, intensity) %>% #put data back into wide format
        group_by(m.z) %>% #t test will be applied for each of these metabolites
        dplyr::mutate(pvalue = t.test(unlist(groupOne), unlist(groupTwo))$p.value)
      
      
      ptable$groupOne <- NULL
      ptable$groupTwo <- NULL
      
      ptable
      
      
    })
    
    output$ptable <- renderTable({
      
      ptable <- req(ptable())
      ptable <- merge(key[,1:2], ptable, by.y = "m.z", by.x ="key")
      ptable <- ptable[order(ptable$pvalue),] #order by pvalue
      ptable$key <- NULL
      names(ptable) <- c("m.z", "pvalue")
      ptable <- ptable[!duplicated(ptable),]
      
      head(ptable)
      
      
      })
    
    
    output$downloadPtable <- downloadHandler(
      filename = function(){paste(input$plotTitles, "- ptable.csv")}, 
      content = function(fname){
        
        ptable <- req(ptable())
        ptable <- merge(key[,1:2], ptable, by.y = "m.z", by.x ="key")
        ptable <- ptable[order(ptable$pvalue),] #order by pvalue
        ptable$key <- NULL
        names(ptable) <- c("m.z", "pvalue")
        ptable <- ptable[!duplicated(ptable),]
        write.csv(ptable, fname, row.names = FALSE)
        
      }
    )
    
    output$heatmap <- renderPlot({
      
      ptable <- req(ptable())
      ptable <- ptable[order(ptable$pvalue),]
      #ptable <- ptable[which(ptable$pvalue < 0.05), 1] %>% sapply(as.character)
      ptable <- ptable[1:input$featureNumber, 1] %>% sapply(as.character)
      
      #create a df with the top 100 metabolites between mean Control and AD
      heatmap_data <- transdf[,c(which(as.character(names(transdf)) %in% ptable == TRUE | as.character(names(transdf)) == "Group"))]
      
      rowAnnot <- data.frame(Group = heatmap_data$Group)
      heatmap_data$Group <- NULL
      rownames(rowAnnot) <- rownames(heatmap_data)
      levels(rowAnnot$Group) <- levels(myData$groupIdentities$Group)
      
      Group <- c("#00425A", "#FC7300")
      names(Group) <- c(levels(myData$groupIdentities$Group)[1], levels(myData$groupIdentities$Group)[2])
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
      
      
    }, height = 400)
    
  })
  
  observeEvent(input$limmaButton, {
    
    top.table <- reactive({
      
      req(transdf)
      groupIdentities <- req(myData$groupIdentities)
      
      #Expression data
      assayData <- as.matrix(t(transdf[, 1:(length(transdf)-1)]))
      
      #phenotype data
      pData <- data.frame(Group = transdf$Group)
      
      
      #all(rownames(pData)==colnames(assayData))
      pData <- new("AnnotatedDataFrame", data=pData)
      
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
      
      #Sanity check
      #colSums(design)
      #table(transdf[,"Group"])
      
      # Fit the model
      fit <- lmFit(eSet, design)
      
      # Calculate the t-stat 
      fit <- eBayes(fit)
      
      # Summarize the results
      #results <- decideTests(fit[,"GroupgroupTwo"])
      
      top.table <- topTable(fit, sort.by = "P", n = Inf)
      top.table <- merge(key, top.table, by = "key")
      top.table$key <- NULL
      top.table <- top.table[order(top.table$P),]

    })
    
    output$toptable <- renderTable({
      
      top.table <- req(top.table())
      head(top.table)
      
    })
    
    
    output$downloadToptable <- downloadHandler(
      filename = function(){paste(input$plotTitles, "- top table.csv")}, 
      content = function(fname){
        
        write.csv(top.table(), fname, row.names = FALSE)
        
      }
    )
    
    output$volcanoPlot <- renderPlot({
      
      top.table <- req(top.table())
      
      top.table$difffex <- "NO"
      top.table$difffex[top.table$logFC > 0.6 & top.table$adj.P.Val < 0.05] <- "UP"
      # if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
      top.table$difffex[top.table$logFC < -0.6 & top.table$adj.P.Val < 0.05] <- "DOWN"
      
      top.table$label[top.table$difffex != "NO"] <- top.table$m.z[top.table$difffex != "NO"]
      
      ggplot(data=top.table, aes(x=logFC, y=-log10(adj.P.Val), color=difffex, label = label)) +
        geom_point() +
        geom_label_repel(max.overlaps=20) +
        theme_bw() +
        scale_color_manual(values=c("blue", "black", "red")) +
        geom_vline(xintercept=c(-0.6, 0.6), col="red") +
        geom_hline(yintercept=-log10(0.05), col="red") + 
        labs(x= "Log2 of Fold Change", 
             y= "Significance (-log10P)",
             color = "Enrichment")
      
      
    })
    
  })
  
  IPS <- reactive({
    
    req(input$file2)
    
    # Read metaboanalyst pathway hits output
    IPS <- read.csv(input$file2$datapath)
    IPS <- IPS[,1:6] # Keep only relevant columns
    
    # Calculate IPS value
    IPS$Weighted_Metabolites <- ((IPS$Hits.sig + 1) ^ 2) + (IPS$Hits.total - IPS$Hits.sig)
    IPS$Numerator <- IPS$Weighted_Metabolites / (IPS$Pathway.total * IPS$Expected)
    IPS$Squared_P_Value <- IPS$FET ^ 2
    IPS$IPS_Value <- IPS$Numerator / IPS$Squared_P_Value + 1
    
    
    # Clean up (Keep only pathway name, IPS, and Map)
    IPS <- IPS[,c(1,10)]
    IPS$log2IPS <- log2(IPS$IPS_Value)
    IPS <- IPS[order(IPS$IPS_Value), ]
    names(IPS) <- c("Pathway", "IPS", "Log2(IPS)")
    
    IPS <- IPS
    
  })
  
  output$IPS <- renderTable(req(IPS()))
  
  output$downloadIPS <- downloadHandler(
    filename = function(){paste("IPS table.csv")}, 
    content = function(fname){
      
      write.csv(IPS(), fname, row.names = FALSE)
      
    }
  )
  
  
}




