server <- function(input, output, session) {
  
  options(shiny.maxRequestSize=100*1024^2)
  myData <- reactiveValues()
  
  
  # Preview Data ----
  
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
  
  # Group identities ---
  
  output$Groups <- renderText({ 
    
    groupIdentities <- req(myData$groupIdentities)
    
    updateSelectInput(session = session, inputId = "relevelSelector", choices = levels(groupIdentities$Group))
      
    paste("<p><b>Control:</b>", levels(groupIdentities$Group)[1], "</p>", 
          "<p><b>Experimental:</b>", 
          paste(levels(groupIdentities$Group)[2:length(levels(groupIdentities$Group))], collapse = ", "), 
          "</p>")
    
  })
  
  
  # Relevel ----
  
  # Need to fix this so that people can relevel after nmds is generated and outliers are removed
  observeEvent(input$relevel, {
    
      req(input$file1)
      groupIdentities <- myData$groupIdentities
      groupIdentities$Group <- relevel(groupIdentities$Group, input$relevelSelector)
      myData$groupIdentities <- groupIdentities
      
    
  })
  
  
  observeEvent(input$Submit, {
    
    showModal(modalDialog("Processing data...", footer=NULL))
    
    data.scores <- isolate({
      
      groupIdentities <- req(myData$groupIdentities)
      
      # Transpose Data and set column names as m/z
      transdfscores <- as.data.frame(t(workingDf))
      
      # Turn all the intensities into numerics and format as numeric matrix
      transdfscores <- sapply(transdfscores, as.numeric)
      colnames(transdfscores) <- workingDf$m.z
      transdfscores <- transdfscores[-1,]
      rownames(transdfscores) <- rownames(groupIdentities)
      
      # Run NMDS
      set.seed(123)
      nmds <- metaMDS(transdfscores, distance = "bray")
      
      # Extract NMDS scores (x and y coordinates)
      data.scores <- as.data.frame(scores(nmds)$sites)
      
    })
    
    isolate({
      
      data.scores <- req(data.scores)
      
      if(input$outlierRemove == TRUE){
        
        # Identify outlier samples based on nmds values and safe a copy of them in directory
        outliers <- rbind(data.scores[which((data.scores$NMDS1 %in% boxplot.stats(data.scores$NMDS1)$out)), ],
                          data.scores[which((data.scores$NMDS2 %in% boxplot.stats(data.scores$NMDS2)$out)), ])
        # write.csv(outliers, paste(i, "- Outliers"), row.names = TRUE)
        
        # remove outliers from all data frames
        data.scores <- data.scores[which(!(rownames(data.scores) %in% rownames(outliers))), ]
        groupIdentities <- data.frame(Group = groupIdentities[which(!(rownames(groupIdentities) %in% rownames(outliers))), ])
        rownames(groupIdentities) <- rownames(data.scores)
        workingDf <- workingDf[, which(!(colnames(workingDf) %in% rownames(outliers)))]
        
        myData$outliers <- rownames(outliers)
        workingDf <<- workingDf
        myData$groupIdentities  <- groupIdentities
        
      }
      
    })
    
    output$outliers <- renderText(
      
      if(input$outlierRemove == TRUE){
        paste(
          "<p>The following samples were idenetified as outliers and have been removed</p> <ul> <li>",
          paste(myData$outliers, collapse='</li><li>'),
          "</li></ul>"
        ) 
      }
    )
    
    
    # P-table ----
    
    ptable <- reactive({
      
      req(input$file1)
      req(workingDf)
      
      groupIdentities <- myData$groupIdentities
      
      
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
      groupIdentities[which(groupIdentities$Group == levels(groupIdentities$Group)[1]),]$AutoGroup <- "group1"
      groupIdentities[which(groupIdentities$Group == levels(groupIdentities$Group)[2]),]$AutoGroup <- "group2"
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
        dplyr::mutate(pvalue = t.test(unlist(group1), unlist(group2))$p.value)
      
      
      ptable$group1 <- NULL
      ptable$group2 <- NULL
      
      ptable
      
      
    })
    
    
    # Preview P-table
    
    output$ptable <- renderTable({
      
      req(input$file1)
      ptable <- req(ptable())
      ptable <- merge(key[,1:2], ptable, by.y = "m.z", by.x ="key")
      ptable <- ptable[order(ptable$pvalue),] #order by pvalue
      ptable$key <- NULL
      names(ptable) <- c("m.z", "pvalue")
      ptable <- ptable[!duplicated(ptable),]
      
      head(ptable)
      
      
    }, digits = -1)
    
    
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
    
    
    # Print NMDS ----
    
    if ("NMDS" %in% input$Plots) {
    
      # input$ANOSIM is not passing correctly into the nmds module
      anoVal <- isolate({ if (input$ANOSIM == TRUE) { anoVal <- TRUE } else { anoVal <- FALSE } })
      generateNmdsServer("nmdsMod", data.scores, myData$groupIdentities, transdf, anoVal, input$plotTitles)
    
    } else { hideTab("plots", target = "NMDS") }
    
    
    
    # Print Heatmap ----
    
    if ("heatmap" %in% input$Plots) {
      
      generateHeatmapServer("heatmapMod", myData$groupIdentities, transdf, ptable())
      
    } else { hideTab("plots", target = "Heatmap") }
    
    
    
    # Run Limma ----
    
    if ("limma" %in% input$Plots) {
      
      generateLimmaServer("limmaMod", transdf, input$plotTitles)
      
    } else { hideTab("plots", target = "Limma") }
    
    
    
    # Run Lasso ----
    
    if ("lasso" %in% input$Plots) {
      
      generateLassoServer("lassoMod", transdf, input$plotTitles)
      
    } else { hideTab("plots", target = "Lasso") }
    
    
    removeModal()
    
  })
  
  
  # IPS outputs -----------------------
  
  IPS <- reactive({
    
    req(input$files2)
    IPS <- data.frame(Pathway = character(),
                      IPS = numeric(),
                      `Log2(IPS)` = numeric(),
                      ID = factor())
    
    
    for(i in 1:length(input$files2[,1])) {
      
      IPSworking <- read.csv(input$files2[i,]$datapath)
      IPSworking <- IPSworking[,1:6] # Keep only relevant columns
      IPSworking$Weighted_Metabolites <- ((IPSworking$Hits.sig + 1) ^ 2) + (IPSworking$Hits.total - IPSworking$Hits.sig)
      IPSworking$Numerator <- IPSworking$Weighted_Metabolites / (IPSworking$Pathway.total * IPSworking$Expected)
      IPSworking$Squared_P_Value <- IPSworking$FET ^ 2
      IPSworking$IPS_Value <- IPSworking$Numerator / IPSworking$Squared_P_Value + 1
      
      # Clean up (Keep only pathway name, IPS, and Map)
      IPSworking <- IPSworking[,c(1,10)]
      IPSworking$log2 <- log2(IPSworking$IPS_Value)
      names(IPSworking) <- c("Pathway", "IPS", "Log2(IPS)")
      IPSworking$ID <- input$files2[i,]$name
      
      IPS <- rbind(IPS, IPSworking)
      
    }
    
    
    IPS <- IPS
    
    
  })
  
  output$IPS <- renderTable({
    req(IPS())
    head(IPS())
    
  })
  
  output$downloadIPS <- downloadHandler(
    filename = function(){paste("IPS table.csv")}, 
    content = function(fname){
      
      write.csv(IPS(), fname, row.names = FALSE)
      
    }
  )
  
  
  output$IPSHeatmap <- renderPlot({
    
    req(IPS())
    IPS <- IPS()
    
    IPS$IPS <- NULL
    
    #Long form
    IPS <- IPS %>%
      pivot_wider(names_from = ID,
                  values_from = `Log2(IPS)`) %>%
      as.data.frame()
    
    rownames(IPS) <- IPS$Pathway
    IPS$Pathway <- NULL
    
    
    pheatmap(IPS,
             cluster_rows = F,
             cluster_cols = F,
             angle_col = 45,
             fontsize_col = 8,
             scale = "column",
    )
    
    
  }, height = 800)
  
  
  
}




