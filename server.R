server <- function(input, output, session) {
  
  options(shiny.maxRequestSize=100*1024^2)
  myData <- reactiveValues()
  
  # Preview Data ----
  
  output$viewInputTable <- renderTable({
    
    req(input$file1)
    workingDf <- read.csv(input$file1$datapath)
    
    if (input$metAnnotations == "annotOnly") {
      
      workingDf <- workingDf[c(1 , which(workingDf$Name != "")),]
      workingDf$m.z <- workingDf$Name
      
    } 
    

    workingDf$m.z <- gsub(" Da", "", workingDf$Bucket.label)
    workingDf$m.z <- as.numeric(workingDf$m.z)
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
  
  
  observeEvent(input$selectGroups, {
    
    groupIdentities <- myData$groupIdentities
    showModal(modalDialog(fluidRow(
      checkboxGroupInput(
        "groupsIncluded",
        "Select groups to analyze:",
        c(sort(unique(groupIdentities$Group))),
      ),
      actionButton("saveGroups", "Save Selection")
    ),
    easyClose = TRUE,
    footer = NULL))
    
  })
  
  observeEvent(input$saveGroups, {

    groupIdentities <- myData$groupIdentities
    groupIdentities <- groupIdentities[which(groupIdentities$Group %in% input$groupsIncluded), , drop = FALSE]
    groupIdentities$Group <- as.factor(as.character(groupIdentities$Group))
    myData$groupIdentities <- groupIdentities
    workingDf <<- workingDf[,c("m.z", rownames(groupIdentities))]
    updateSelectInput(session = session, inputId = "relevelSelector", choices = levels(groupIdentities$Group))
    
  })

  
  observeEvent(input$Submit, {
    showModal(modalDialog("Processing data...", footer=NULL))
    
    data.scores <- isolate({
      browser
      groupIdentities <- req(myData$groupIdentities)
      
      # Transpose Data and set column names as m/z
      transdfscores <- as.data.frame(t(workingDf))
      
      # Turn all the intensities into numerics and format as numeric matrix
      transdfscores <- transdfscores[-1,]
      transdfscores <- sapply(transdfscores, as.numeric)
      colnames(transdfscores) <- workingDf$m.z
      rownames(transdfscores) <- rownames(groupIdentities)
      
      # Run NMDS
      set.seed(123)
      nmds <- metaMDS(transdfscores, distance = "bray")
      
      # Extract NMDS scores (x and y coordinates)
      data.scores <- as.data.frame(scores(nmds)$sites)
      
    })
    
    isolate({
      
      data.scores <- req(data.scores)
      
      if (input$outlierRemove == TRUE) {
        
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
    
    
    transdf <- isolate({
      
      
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
      
      for (i in seq(1:length(levels(groupIdentities$Group)))) {
        
        groupIdentities[which(groupIdentities$Group == levels(groupIdentities$Group)[i]), ]$AutoGroup <- paste0("group", i)
        
      }
      
      transdf$Group <- groupIdentities$AutoGroup
      myData$groupIdentities <- groupIdentities
      transdf$Group <- as.factor(transdf$Group)
      transdf <- transdf[, c(which(colSums(transdf[, 1:(length(transdf)-1)]) != 0), length(transdf))]
      transdf <<- transdf
      
    })
    
    # P-table ----
    
    
    
    ptable <- reactive({
      
      #Create a table of pvalues 
      ptable <- transdf %>% #select only a few columns to start with
        gather(key = m.z, value = intensity, -Group) %>% #melt the data into a long format
        group_by(Group, m.z) %>% #group intensities by metabolite and cohort
        dplyr::summarize(intensity = list(intensity)) %>% #put the intensities for all samples in a cohort for a particular metabolite in one column
        spread(Group, intensity)
      
      
      
      for (i in seq(1:(length(levels(groupIdentities$Group)) - 1))) {
        
        workingptable <- ptable %>% #put data back into wide format
          group_by(m.z) %>% #t test will be applied for each of these metabolites
          dplyr::mutate(
            pvalue = t.test(
                     unlist(get(paste0("group", 1))), 
                     unlist(get(paste0("group", i+1)))
                   )$p.value 
            )

        vaName <- paste0("pvalue - ", levels(groupIdentities$Group)[i+1])
        ptable[[vaName]] <- workingptable$pvalue
        ptable <- ptable %>% dplyr::select(-paste0("group", i+1))
        
      }
      
      ptable$group1 <- NULL
      
      ptable
      
      
    })
    
    
    # Preview P-table
    
    output$ptable <- renderTable({
      
      req(input$file1)
      ptable <- req(ptable())
      ptable <- merge(key[,1:2], ptable, by.y = "m.z", by.x ="key")
      ptable$key <- NULL
      ptable <- ptable[order(ptable[, 2]),] #order by pvalue
      ptable <- ptable[!duplicated(ptable),]
      names(ptable) <- gsub("pvalue - ", "", names(ptable))
      
      head(ptable)
      
      
    }, digits = -1)
    
    temp_directory <- paste0(tempdir(), as.integer(Sys.time()), "/")
    
    output$downloadPtable <- downloadHandler(
      filename = function() {
        paste0(input$plotTitles, "- ptables.zip")
      },
      
      content = function(file) {
        ptable = as.data.frame(req(ptable()))
        ptable = merge(key[,1:2], ptable, by.y = "m.z", by.x ="key")
        ptable$key = NULL
        
        dir.create(temp_directory)
        
        for (i in names(ptable)) {
          
          if (i == "m.z") {} else {
            p = ptable[, c("m.z", i)]
            names(p) = c("m.z", "pvalue")
            p = p[which(is.na(p[, 2]) == FALSE),]
            p = p[order(p[, 2]),]
            names(p) = c("m.z", "pvalue")
            file_name = gsub("pvalue - ", "", i)
            write.csv(p, paste0(temp_directory, file_name, ".csv"), row.names = FALSE)
          }
        }
        
        zip::zip(
          zipfile = file,
          files = dir(temp_directory),
          root = temp_directory
        )
        
      },
      contentType = "application/zip"
      
    )

    
    # Print NMDS ----
    
    if ("NMDS" %in% input$Plots) {
      
      req(transdf)
      # input$ANOSIM is not passing correctly into the nmds module
      anoVal <- isolate({ if (input$ANOSIM == TRUE) { anoVal <- TRUE } else { anoVal <- FALSE } })
      generateNmdsServer("nmdsMod", data.scores, myData$groupIdentities, transdf, anoVal, input$plotTitles)
    
    } else { hideTab("plots", target = "NMDS") }
    
    
    
    # Print Heatmap ----
    
    if ("heatmap" %in% input$Plots) {
      
      req(transdf)
      generateHeatmapServer("heatmapMod", myData$groupIdentities, transdf, ptable())
      
    } else { hideTab("plots", target = "Heatmap") }
    
    
    
    # Run Limma ----
    
    if ("limma" %in% input$Plots) {
      
      req(transdf)
      generateLimmaServer("limmaMod", transdf, input$plotTitles, myData$groupIdentities)
      
    } else { hideTab("plots", target = "Limma") }
    
    
    
    # Run Lasso ----
    
    if ("lasso" %in% input$Plots) {
      
      req(transdf)
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
      IPSworking$ID <- gsub(".csv", "", input$files2[i,]$name)
      IPSworking <- IPSworking[order(-IPSworking$IPS),]
      
      IPS <- rbind(IPS, IPSworking)
      
    }
    
    IPS$ID <- as.factor(IPS$ID)
    updateSelectInput(session = session, inputId = "orderBySelector", choices = levels(IPS$ID))
    
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
    
    IPS <- req(IPS())
    
    IPS$IPS <- NULL
    
    #Long form
    IPS <- IPS %>%
      pivot_wider(names_from = ID,
                  values_from = `Log2(IPS)`) %>%
      as.data.frame()
    
    rownames(IPS) <- IPS$Pathway
    IPS$Pathway <- NULL
    
    if (input$naToZero == TRUE) {
      
      IPS[is.na(IPS)] <- 0
      
    }
    
    if (input$orderBySelector != "") {
        
      col <- input$orderBySelector
      IPS <- IPS[order(-IPS[,col]),]
      IPS <- IPS %>%
        dplyr::select(all_of(col), everything())
      
      
      pheatmap(IPS,
               cluster_rows = input$IPSclustRow,
               cluster_cols = input$IPSclustCol,
               angle_col = 45,
               fontsize = 12,
               scale = "column",
      )
      
    }
    
  })
  
  
  
}




