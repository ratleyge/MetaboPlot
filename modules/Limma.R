generateLimmaUI <- function (id) {
  ns <- NS(id)
  
  tagList(h3("Limma & Volcano Plot"),
          uiOutput(ns("LimmaUIOutput")),)
  
}


generateLimmaServer <-
  function(id, transdf, plotTitle, groupIdentities, transformQC, scalingQC) {
    moduleServer(id, function(input, output, session) {
      makeTopTable <- function(transdf) {
        
        #Expression data
        assayData <-
          as.matrix(t(transdf[, 1:(length(names(transdf)) - 1)]))
        
        # Quality control functions
        if (!transformQC) assayData <- logTransform(assayData)
        if (!scalingQC) assayData <- paretoScale(assayData)
        
        #phenotype data
        pData <- data.frame(Group = transdf$Group)
        pData <- new("AnnotatedDataFrame", data = pData)
        sampleNames(pData) <- colnames(assayData)
        
        #Feature Data
        FeatureData <-
          as.data.frame(names(transdf[, 1:(length(transdf) - 1)]))
        names(FeatureData) <- "key"
        FeatureData <- new("AnnotatedDataFrame", data = FeatureData)
        featureNames(FeatureData) <- rownames(assayData)
        
        
        eSet <- ExpressionSet(
          assayData = assayData,
          phenoData = pData,
          featureData = FeatureData
        )
        

        #Create a model matrix
        design <- model.matrix( ~ Group, data = pData(eSet))
        
        # Fit the model
        fit <- lmFit(eSet, design)
        
        # Calculate the t-stat
        fit <- eBayes(fit)
        
        
        # Summarize the results
        #results <- decideTests(fit[,"Groupgroup2"])
        
        top.table <- topTable(fit, sort.by = "none", n = Inf)
        top.table <- merge(key, top.table, by = "key")
        top.table$key <- NULL
        top.table <- top.table[order(top.table$P), ]
        return(top.table)
        
      }
      
      # Generate Volcano ----
      
      output$LimmaUIOutput <- renderUI({
        FileNumber <- length(levels(transdf$Group))
        plot_output_list <- list()
        
        for (i in 2:FileNumber) {
          ns <- session$ns
          
          plot_output_list <- append(plot_output_list, list(htmlOutput(ns(paste0("LimmaTitle", i)))))
          plot_output_list <- append(plot_output_list, list(tableOutput(ns(paste0("toptable", i)))))
          plot_output_list <- append(plot_output_list, list(downloadButton(ns(paste0('downloadToptable', i)), "Download Table")))
          plot_output_list <- append(plot_output_list, list(
            withSpinner(
              plotOutput(ns(paste0(
                "eachVolcano", i
              )), width = "100%", height = "400px"),
              type = getOption("spinner.type", default = 5)
            )
          ))
          
        }
        
        do.call(tagList, plot_output_list)
        
      })
      
      
      
      observe({
        FileNumber <- length(levels(transdf$Group))
        
        for (i in 2:FileNumber) {
          local({
            my_i <- i
            tableName <- paste0("toptable", my_i)
            limmaTitle <- paste0("LimmaTitle", my_i)
            downloadTopName <- paste0('downloadToptable', my_i)
            plotname <- paste0("eachVolcano", my_i)
            
            workingTransDf <-
              transdf[which(
                transdf$Group == levels(transdf$Group)[1] |
                  transdf$Group == levels(transdf$Group)[my_i]
              ),]
            workingTransDf <-
              workingTransDf[, c(which(colSums(workingTransDf[, 1:(length(workingTransDf) -
                                                                     1)]) != 0), length(workingTransDf))]
            workingTransDf$Group <-
              as.factor(as.character(workingTransDf$Group))
            workingTopTable <- makeTopTable(workingTransDf)
            
            output[[limmaTitle]] <-
              renderText(paste0("<h4>", levels(groupIdentities$Group)[my_i], "</h4>"))
            output[[tableName]] <- renderTable(head(workingTopTable))
            
            output[[downloadTopName]] <- downloadHandler(
              filename = function() {
                paste(plotTitle,
                      "-",
                      levels(groupIdentities$Group)[my_i],
                      "- top table.csv")
              },
              content = function(fname) {
                write.csv(workingTopTable, fname, row.names = FALSE)
              }
            )
            
            output[[plotname]] <- renderPlot({
              
              workingTopTable$diffex <- "NO"
              workingTopTable$label <- NA
                
              if (length(rownames(workingTopTable[which(workingTopTable$logFC > 0.6), ])) > 0 & length(rownames(workingTopTable[which(workingTopTable[which(workingTopTable$logFC > 0.6), ]$adj.P.Val < 0.05), ])) > 0) {
                
                workingTopTable[which(workingTopTable$logFC > 0.6 &
                                        workingTopTable$adj.P.Val < 0.05),]$diffex <- "UP"
                
                workingTopTable$label[workingTopTable$diffex != "NO"] <-
                  workingTopTable$m.z[workingTopTable$diffex != "NO"]
                
              }
              
              if (length(rownames(workingTopTable[which(workingTopTable$logFC > 0.6), ])) > 0 & length(rownames(workingTopTable[which(workingTopTable[which(workingTopTable$logFC < -0.6), ]$adj.P.Val < 0.05), ])) > 0) {
                # if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
                workingTopTable[which(workingTopTable$logFC < -0.6 &
                                        workingTopTable$adj.P.Val < 0.05),]$diffex <- "DOWN"
                
                workingTopTable$label[workingTopTable$diffex != "NO"] <-
                  workingTopTable$m.z[workingTopTable$diffex != "NO"]
              }
              

              p <-
                ggplot(data = workingTopTable,
                       aes(
                         x = logFC,
                         y = -log10(adj.P.Val),
                         color = diffex,
                         label = label
                       )) +
                geom_point() +
                geom_label_repel(max.overlaps = 20) +
                theme_bw() +
                scale_color_manual(values = c("DOWN" = "blue", "NO" = "black", "UP" = "red")) +
                geom_vline(xintercept = c(-0.6, 0.6), col = "red") +
                geom_hline(yintercept = -log10(0.05),
                           col = "red") +
                labs(x = "Log2 of Fold Change",
                     y = "Significance (-log10P)",
                     color = "Enrichment")
              
              print(p)
              
            })
            
          })
        }
      })
    })
  }
