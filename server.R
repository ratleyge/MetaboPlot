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
      transDf <- as.data.frame(t(workingDf))
      
      # Turn all the intensities into numerics and format as numeric matrix
      transDf <- sapply(transDf, as.numeric)
      colnames(transDf) <- workingDf$m.z
      transDf <- transDf[-1,]
      rownames(transDf) <- rownames(groupIdentities)
      
      # Run NMDS
      set.seed(123)
      nmds <- metaMDS(transDf, distance = "bray")
      
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
        transDf <- transDf[which(!(rownames(transDf) %in% rownames(outliers))), ]
        workingDf <- workingDf[, which(!(colnames(workingDf) %in% rownames(outliers)))]
        
        myData$outliers <- rownames(outliers)
        workingDf <<- workingDf
        myData$groupIdentities  <- groupIdentities
        
      }
      
    })
    
    # input$ANOSIM is not passing correctly into the nmds module
    anoVal <- isolate({ if (input$ANOSIM == TRUE) { anoVal <- TRUE } else { anoVal <- FALSE } })
    generateNmdsServer("nmdsMod", data.scores, myData$groupIdentities, transDf, anoVal)
    
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
    
    
    
    # Print Heatmap ----
    
    output$heatmap <- renderPlot({
      
      if ("heatmap" %in% input$Plots) {
        
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
      } else {
        
        hideTab("plots", target = "Heatmap")
        
      }
      
      
    })
    
    
    # Run Limma ----
    
    top.table <- isolate({
      
      if ("limma" %in% input$Plots) {
        
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
        #results <- decideTests(fit[,"Groupgroup2"])
        
        top.table <- topTable(fit, sort.by = "P", n = Inf)
        top.table <- merge(key, top.table, by = "key")
        top.table$key <- NULL
        top.table <- top.table[order(top.table$P),]
      } else {
        
        hideTab("plots", target = "Limma")
        
      }
    })
    
    output$toptable <- renderTable({
      if ("limma" %in% input$Plots) {
        
        req(top.table)
        head(top.table)
        
      }
      
    })
    
    
    output$downloadToptable <- downloadHandler(
      filename = function(){paste(input$plotTitles, "- top table.csv")}, 
      content = function(fname){
        
        write.csv(top.table, fname, row.names = FALSE)
        
      }
    )
    
    # Generate Volcano ----
    
    output$volcanoPlot <- renderPlot({
      
      if ("limma" %in% input$Plots) {
        
        req(top.table)
        
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
        
      }
    })
    
    
    
    
    
    # Run Lasso ----
    
    coefs <- isolate({
      
      if ("lasso" %in% input$Plots) {
        
        y = transdf$Group
        x = as.matrix(transdf[, 1:(length(transdf)-1)])
        
        cv_model <- cv.glmnet(scale(x), 
                              y, 
                              alpha = 1, 
                              standardize = TRUE,
                              nfolds = 10,
                              family = "binomial")
        
        #find optimal lambda value that minimizes test MSE
        best_lambda <- cv_model$lambda.min
        
        #produce plot of test MSE by lambda value
        #plot(cv_model)
        #plot(cv_model$glmnet.fit, col=blues9, xvar="lambda", label=F)
        #abline(v=log(best_lambda), col="red")
        #abline(v=log(cv_model$lambda.1se), col="red")
        
        best_model <- glmnet(scale(x), 
                             y, 
                             alpha = 1, 
                             standardize = TRUE,
                             lambda = best_lambda, 
                             family = "binomial")
        
        
        coef <- coef(best_model, s = best_model$lambda.min)
        coefs <- data.frame("key" = coef@Dimnames[[1]][2:length(coef)], "Beta" = coef[2:length(coef)])
        coefs <- merge(key, coefs, by = "key")
        coefs$key <- NULL
        coefs <- coefs[order(-coefs$Beta),]
        
      }
      
      else {
        
        hideTab("plots", target = "Lasso")
        
      }
      
    })
    
    output$LassoTable <- renderTable({
      
      if ("lasso" %in% input$Plots) {
        
        req(coefs)
        head(coefs)
        
      }
      
    }, digits = -1)
    
    output$downloadLasso <- downloadHandler(
      filename = function(){paste(input$plotTitles, "- Lasso Coefficients.csv")}, 
      content = function(fname){
        req(coefs)
        write.csv(coefs, fname, row.names = FALSE)
        
      }
    )
    
    output$LassoCoefficients <- renderPlot({
      
      
      if ("lasso" %in% input$Plots) {
        
        req(coefs)
        
        coefs <- coefs[abs(coefs$Beta) > 0, ]
        coefs <- coefs[order(-abs(coefs$Beta)), ]
        
        g <- (ggplot(coefs, aes(x = reorder(m.z, Beta), y = Beta)) +
                geom_bar(stat = "identity", 
                         position = position_stack(),
                         color = "white", 
                         fill = "lightblue") +
                ggtitle(input$plotTitles) +
                coord_flip() +
                ylab("Beta Coefficients") +
                xlab("Input Variables") +
                theme(plot.margin = unit(c(0.5,1.5,0.5,0.5), "cm")) +
                theme(text = element_text(size = 15))) +
          geom_hline(yintercept = 0)
        
        g
        
      }
      
    })
    
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




