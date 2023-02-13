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
  
  observeEvent(input$relevel, {
    
    groupIdentities <- req(myData$groupIdentities)
    groupIdentities$Group <- relevel(groupIdentities$Group, levels(groupIdentities$Group)[2])
    myData$groupIdentities <- groupIdentities
    
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
      
      ptable <- merge(key[,1:2], ptable, by.y = "m.z", by.x ="key")
      ptable <- ptable[order(ptable$pvalue),] #order by pvalue
      ptable$key <- NULL
      names(ptable) <- c("m.z", "pvalue")
      ptable <- ptable[!duplicated(ptable),]
      
      
    })
    
    output$ptable <- renderTable(head(req(ptable())))
    output$downloadPtable <- downloadHandler(
      filename = function(){paste(input$plotTitles, "- ptable.csv")}, 
      content = function(fname){
        write.csv(ptable(), fname, row.names = FALSE)
      }
    )
    
  })
}




