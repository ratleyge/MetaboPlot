server <- function(input, output, session) {
  
  options(shiny.maxRequestSize=100*1024^2)
  
  myData <- reactiveValues()
  
  # Outputs 
  output$viewInputTable <- renderTable({
    
    req(input$file1)
    df <- read.csv(input$file1$datapath)
    
    df$Bucket.label <- NULL
    df$RT <- NULL
    df$Name <- NULL
    df$Formula <- NULL
    myData$df <- df
    
    groupIdentities <- data.frame(df[1, ])
    groupIdentities$m.z <- NULL
    groupIdentities <- as.data.frame(t(groupIdentities))
    names(groupIdentities) <- "Group"
    groupIdentities$Group <- as.factor(groupIdentities$Group)
    myData$groupIdentities <- groupIdentities
    
    head(df[-1, 1:6])
    
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
      
      df <- req(myData$df)
      groupIdentities <- req(myData$groupIdentities)
      
      showModal(modalDialog("Generating NMDS plot. This may take a moment.", footer=NULL))
      
      df <- df[-1,]
      
      # Transpose Data and set column names as m/z
      transdf <- as.data.frame(t(df))
      
      # Turn all the intensities into numerics and format as numeric matrix
      transdf <- sapply(transdf, as.numeric)
      colnames(transdf) <- df$m.z
      transdf <- transdf[-1,]
      rownames(transdf) <- rownames(groupIdentities)
      
      # Run NMDS
      set.seed(123)
      nmds <- metaMDS(transdf, distance = "bray")
      nmds
      
      # Extract NMDS scores (x and y coordinates)
      data.scores <- as.data.frame(scores(nmds)$sites)
      
      # Identify outlier samples based on nmds values and safe a copy of them in directory
      outliers <- rbind(data.scores[which((data.scores$NMDS1 %in% boxplot.stats(data.scores$NMDS1)$out)), ],
                        data.scores[which((data.scores$NMDS2 %in% boxplot.stats(data.scores$NMDS2)$out)), ])
      # write.csv(outliers, paste(i, "- Outliers"), row.names = TRUE)
      
      # remove outliers from all data frames
      data.scores <- data.scores[which(!(rownames(data.scores) %in% rownames(outliers))), ]
      groupIdentities <- data.frame(Group = groupIdentities[which(!(rownames(groupIdentities) %in% rownames(outliers))), ])
      rownames(groupIdentities) <- rownames(data.scores)
      transdf <- transdf[which(!(rownames(transdf) %in% rownames(outliers))), ]
      df <- df[, which(!(colnames(df) %in% rownames(outliers)))]
      
      
      
      
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
              legend.key=element_blank()) + 
        labs(x = "NMDS1", colour = "Group", y = "NMDS2", shape = "Group", title = paste(input$plotTitles, "NMDS")) 
      
      
      if (input$ANOSIM == TRUE) {
        
        ano <- anosim(transdf, groupIdentities[, 1], distance = "bray", permutations = 9999)
        xx <- xx + labs(subtitle = paste("ANOSIM stat:", 
                                         signif(ano$statistic, digits = 2),
                                         "Significance:", 
                                         signif(ano$signif, digits = 2)))
        
      }
      
      
      
      myData$outliers <- rownames(outliers)
      myData$df_OutliersRemoved <- df
      myData$groupIdentities_OutliersRemoved  <- groupIdentities
      myData$transdf_OutliersRemoved  <- transdf
      
      removeModal()
      print(xx)
      
    })
    
    
    output$outliers <- renderText(
      paste(
        "The following samples were idenetified as outliers and are not pictured in the NMDS plot:",
        paste(myData$outliers, collapse=', '),
        "These samples will be removed from further analysis."
      )
    )
    
  })
}