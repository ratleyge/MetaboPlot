
createNMDS <- function(id, groupIdentities, ) {
  
    moduleServer(
      
      id,
      
      ## Below is the module function
      function(input, output, session) {
        
        # The selected file, if any
        userFile <- reactive({
          # If no file is selected, don't do anything
          validate(need(input$file, message = FALSE))
          input$file
        })
        
        # The user's data, parsed into a data frame
        dataframe <- reactive({
          read.csv(userFile()$datapath,
                   header = input$heading,
                   quote = input$quote,
                   stringsAsFactors = stringsAsFactors)
        })
        
        # We can run observers in here if we want to
        observe({
          msg <- sprintf("File %s was uploaded", userFile()$name)
          cat(msg, "\n")
        })
        
        # Return the reactive that yields the data frame
        return(dataframe)
      }
    )    
  }




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
  
  if(input$outlierRemove == TRUE) {
    
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
  
  
  # Create NMDS plot
  
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
