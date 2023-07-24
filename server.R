server <- function(input, output, session) {

  options(shiny.maxRequestSize=200*1024^2) # Allow users to upload large files
  myData <- reactiveValues() # This allows you to define variables that will carry across functions
  
  # Preview Data ----

    # This function reads in the intensity table which the user has uploaded
    # It also creates several variables including: workingDf, groupIdentities, and key
    output$viewInputTable <- renderTable({
      
      req(input$file1)
      workingDf <- read.csv(input$file1$datapath) # Read the file
      
      # remove rows with na's, MeatboScape will sometimes have empty rows with all na at the bottom
      workingDf <- rbind(workingDf[1, ], na.omit(workingDf)) 
      
      if (input$metAnnotations == "annotOnly") {
        
        # If the user selects annotations only, then the value of m/z should come from the name column
        workingDf <- workingDf[c(1 , which(workingDf$Name != "")),]
        workingDf$m.z <- workingDf$Name
        
      } else if (input$metAnnotations == "mzOnly") {
        
        # If the user selects m/z only then the value of mz should come from the bucket label column
        workingDf$m.z <- gsub(" Da", "", workingDf$Bucket.label)
        workingDf$m.z <- as.numeric(workingDf$m.z)
        
      } else if (input$metAnnotations == "mzAnnot") {
        
        workingDf[which(!is.na(workingDf$Name)), ]$m.z <- workingDf[which(!is.na(workingDf$Name)), ]$Name
        workingDf[which(is.na(workingDf$Name)), ]$m.z <- gsub(" Da", "", workingDf[which(is.na(workingDf$Name)), ]$Bucket.label)
        
      }
      
      # Get rid of useless columns
      workingDf$Bucket.label <- NULL
      workingDf$RT <- NULL
      workingDf$Name <- NULL
      workingDf$Formula <- NULL
      
      # Create a data frame that matches group identity in row 1 to each sample
      groupIdentities <- data.frame(workingDf[1, ])
      groupIdentities$m.z <- NULL
      groupIdentities <- as.data.frame(t(groupIdentities))
      names(groupIdentities) <- "Group"
      groupIdentities$Group <- as.factor(groupIdentities$Group)
      workingDf <- workingDf[-1,] # remove the group row from working df when finished
      
      # Sometimes multiple rows have the same name, create a dummy variable and store the names in the key data frame
      rownames(workingDf) <- paste("met", 1:length(workingDf$m.z))
      key <- data.frame(Key = paste("met", 1:length(workingDf$m.z)), m.z = workingDf$m.z)
      workingDf$m.z <- NULL 
      
      # The working Df should only contain intensity values, sample names will be in the column names, key for the metabolite should be in the row names
      workingDf[] <- sapply(workingDf, as.numeric)

      # Remove rows with a certain % missingness as defined by the user
      workingDf <- removeGroupMissingness(workingDf, groupIdentities, input$wholeMissingness, input$groupMissingness)
      
      # If a sample has all 0s, remove it
      workingDf <- workingDf[, which(colSums(workingDf) > 0), drop = FALSE]
      groupIdentities <- groupIdentities[names(workingDf),, drop = FALSE] # be consistent with groupIdenetities
      
      # Quality control functions
      if (input$transformQC) workingDf <- logTransform(workingDf)
      if (input$scalingQC) workingDf <- paretoScale(workingDf)
      
      # Reset the row names to be the values of m/z for display purposes
      workingDf$Key <- rownames(workingDf)
      workingDf <- merge(key, workingDf, by = "Key")
      workingDf$Key <- NULL
      
      # Make working Df and key accessible in other functions
      workingDf <<- workingDf
      myData$groupIdentities <- groupIdentities
      
      # Display the first six rows
      head(workingDf[, 1:6])
      
    }) 
  
  output$dimensionOutput <- renderText({
    
    # Any time one of these inputs changes, recalculate the data dimensions
    req(input$file1)
    req(input$metAnnotations)
    req(input$wholeMissingness)
    req(input$groupMissingness)
    req(myData$groupIdentities)
    
    paste0(length(rownames(workingDf)), " peaks across ", length(names(workingDf)), " samples.")
    
  })
    
  
  # Group identities ----
    
    # This will display the group identities in the tab panel
    output$Groups <- renderText({ 
      
      groupIdentities <- req(myData$groupIdentities)
      updateSelectInput(session = session, inputId = "relevelSelector", choices = levels(groupIdentities$Group))
      paste("<p><b>Control:</b>", levels(groupIdentities$Group)[1], "</p>", 
            "<p><b>Experimental:</b>", 
            paste(levels(groupIdentities$Group)[2:length(levels(groupIdentities$Group))], collapse = ", "), 
            "</p>")
      
    })
  
  
  # Relevel and Select Groups ----
  
    # Need to fix this so that people can relevel after nmds is generated and outliers are removed
    # This allows the user to select a control group for later analyses
    observeEvent(input$relevel, {
      
      req(input$file1)
      groupIdentities <- myData$groupIdentities
      groupIdentities$Group <- relevel(groupIdentities$Group, input$relevelSelector)
      myData$groupIdentities <- groupIdentities
      
    })
    
    # When a user selects the button "Select Groups" 
    # it allows them to subset their data based on group identity in a popup window
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
    
    # Once the user clicks the save groups button, it removes the groups from the workingdf
    observeEvent(input$saveGroups, {
      
      groupIdentities <- myData$groupIdentities
      groupIdentities <- groupIdentities[which(groupIdentities$Group %in% input$groupsIncluded), , drop = FALSE]
      groupIdentities$Group <- as.factor(as.character(groupIdentities$Group))
      myData$groupIdentities <- groupIdentities
      workingDf <<- workingDf[,c("m.z", rownames(groupIdentities))]
      updateSelectInput(session = session, inputId = "relevelSelector", choices = levels(groupIdentities$Group))
      
    })
    
  
  # Submit Data ----
    
    # When the user submits their data, the code will make ptables and perform selected analyses
    observeEvent(input$Submit, {
      
      # Dealing with outliers using nmds ----
        
        # Show a modal while analyses are happening
        # I should probably fix this because it always disappears too early
        showModal(modalDialog("Processing data...", footer=NULL))
        
        # Data scores runs nmds (even if the user doesn't select it)
        # so that you can remove outliers based on the results
        data.scores <- isolate({
          
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
        
        # If the user chase to remove outliers, this will remove them based on nmds values
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
        
        # Display which samples were removed
        output$outliers <- renderText(
          
          if(input$outlierRemove == TRUE){
            paste(
              "<p>The following samples were idenetified as outliers and have been removed</p> <ul> <li>",
              paste(myData$outliers, collapse='</li><li>'),
              "</li></ul>"
              
            ) 
          }
        )
       
      # Transpose data frame  ----   
        
        # Create the variable transdf in which m/z values are column names and samples are rows
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
          
          
          if (length(levels(groupIdentities$Group)) > 2) {
            
            AnovaTable <- CalculateAnova(transdf)
            ptable <- merge(ptable, AnovaTable, by = "m.z")
            
          }
          
          
          ptable
          
          
        })
        
        
        # Preview P-table
        
        output$ptable <- renderTable({
          
          req(input$file1)
          ptable <- req(ptable())
          ptable <- merge(key[,1:2], ptable, by.y = "m.z", by.x ="key")
          ptable$key <- NULL
          
          if (length(levels(groupIdentities$Group)) > 2) {
            
            ptable <- ptable[order(ptable$AnovaPvalue),]
            
          } else {
            
            ptable <- ptable[order(ptable[, 2]),] #order by pvalue
            
          }
          
          
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
        generateLimmaServer("limmaMod", transdf, input$plotTitles, myData$groupIdentities, input$transformQC, input$scalingQC)
        
      } else { hideTab("plots", target = "Limma") }
      
      
      
      # Run Lasso ----
      
      if ("lasso" %in% input$Plots) {
        
        req(transdf)
        generateLassoServer("lassoMod", transdf, input$plotTitles, myData$groupIdentities)
        
      } else { hideTab("plots", target = "Lasso") }
      
      
      removeModal()
      
  })
    
    # MetaboAnalyst outputs ---------------------
    
    output$metaboAnalystInstalled <- renderText({
      if(!require(MetaboAnalystR)) {
        paste("<h1>You do not have MetaboAnalystR installed on your computer, so you will not be able to use the functions in this tab.</h1>")
      }
    })
    
    # Preview Data ----
    #isolate function- runs one time after hit submit
    observeEvent(input$SubmitMA, 
                 { 
                   #browser()
                   show_modal_spinner()
                   tempdirect <<- paste0(tempdir(), "FunctionalAnalysis", "/")
                   unlink(paste0(tempdirect), recursive = TRUE)
                   isolate({
                     
                     trace("download.file", where = asNamespace("MetaboAnalystR"), tracer = quote(method <- "wininet"))
                     myOrganism <- input$Organisms
                     myMode <- input$IonMode
                     tempdirect <- paste0(tempdir(), "FunctionalAnalysis", "/")
                     dir.create(tempdirect)
                     print(list.files(tempdirect))
                     plotdataframes <<- data.frame()
                     
                     for (i in 1:length(input$file3$datapath)) {
                       
                       mSet <- NULL
                       #browser()
                       all.mzsn <- NULL
                       mdata.all <- NULL
                       mdata.siggenes <- NULL
                       anal.type <- NULL
                       api.base <- NULL
                       err.vec <- NULL
                       meta.selected <- NULL
                       module.count <- NULL
                       msg.vec <- NULL
                       smpdbpw.count <- NULL
                       url.pre <- NULL
                       
                       
                       getwd()
                       wd_mum_raw<-paste0(getwd(), "/mum_raw.qs")
                       wd_mum_res<-paste0(getwd(), "/mum_res.qs")
                       #browser()
                       wd_hsa_mfn<-paste0(getwd(), "/", input$Organisms,".qs")
                       wd_mummi_matched<-paste0(getwd(), "/mummichog_matched_compound_all.csv")
                       wd_mummi_pathway<-paste0(getwd(), "/mummichog_pathway_enrichment.csv")
                       wd_mummi_query<-paste0(getwd(), "/mummichog_query.json")
                       wd_peaks<-paste0(getwd(), "/peaks_to_paths_0_dpi72.png")
                       wd_scattermum<-paste0(getwd(), "/scattermum.json")
                       
                       
                       #Create dynamic file name and write output as csv
                       
                       subfolder_names<- input$file3[i,]$name
                       directoryName <<- paste0(tempdirect, gsub("_pval.csv","",input$file3[i,]$name), "/")
                       dir.create(directoryName)
                       mSet<-InitDataObjects("mass_all", "mummichog", FALSE);
                       mSet<-SetPeakFormat(mSet, "rmp")
                       mSet<-UpdateInstrumentParameters(mSet, input$MT, myMode, "no", 0.02);
                       mSet<-Read.PeakListData(mSet, input$file3[i,]$datapath);
                       mSet<-SanityCheckMummichogData(mSet)
                       mSet<-SetPeakEnrichMethod(mSet, "mum", "v2")
                       mSet<-SetMummichogPval(mSet, 0.05)
                       mSet<-PerformPSEA(mSet, myOrganism, "current", 3 , 100)
                       mSet<-PlotPeaks2Paths(mSet, "peaks_to_paths_0_", "png", 72, width=NA)
                       enrichmentData <- fromJSON(file = wd_scattermum)
                       enrichmentData <- data.frame(
                         Pathway = enrichmentData$pathnames, 
                         Pvalue = enrichmentData$pval, 
                         log10Pvalue = -log10(enrichmentData$pval),
                         Enrichment = enrichmentData$enr
                       )
                       df<-enrichmentData[rev(order(enrichmentData$log10Pvalue)),]
                       df$Label <- NA
                       df[1:min(5,length(rownames(df))),]$Label <- df[1:min(5,length(rownames(df))),]$Pathway
                       ggplot(df, aes(x = Enrichment, y = log10Pvalue, label = Label)) + 
                         geom_point(aes(color = log10Pvalue)) + 
                         scale_color_gradient(low = "Yellow", high = "Red") + 
                         ggrepel::geom_text_repel() + 
                         theme_classic()
                       ggsave(paste0(gsub("_pval.csv","",input$file3[i,]$name), "_myplot.png.png"), path= directoryName)
                       df$sample <- input$file3[i,]$name
                       plotdataframes <<- rbind(plotdataframes, df)
                       
                       
                       
                       
                       file.move(c(wd_mum_raw, wd_mum_res, wd_hsa_mfn,wd_mummi_matched,wd_mummi_pathway,wd_mummi_query, wd_peaks, wd_scattermum), directoryName)
                       print(list.files(tempdirect))
                       #file.remove(from=c(wd_mum_raw, wd_mum_res, wd_hsa_mfn,wd_mummi_matched,wd_mummi_pathway,wd_mummi_query, wd_peaks, wd_scattermum))
                       
                       
                     }
                     
                     
                   })
                   output$downloadfiles <- downloadHandler(
                     filename = function() {
                       paste0(tempdirect, "FunctionalAnalysisResults.zip")
                     },
                     
                     content = function(FunctionalAnalysis2) {
                       folderstodownload<-list.dirs(paste0(tempdirect, gsub("_pval.csv","",input$file3$name)))
                       zip::zipr(
                         zipfile = FunctionalAnalysis2,
                         files= folderstodownload,
                         recurse = TRUE,
                         include_directories = TRUE,
                         root = tempdirect
                       )
                       
                       
                     },
                     contentType = "application/zip"
                     
                   )
                   
                   
                   IPS<- reactive({
                     req(input$file3)
                     IPSinput<-(paste0(tempdirect,gsub("_pval.csv","",input$file3$name), "/mummichog_pathway_enrichment.csv"))
                     IPS <- data.frame(Pathway = character(),
                                       IPS = numeric(),
                                       `Log2(IPS)` = numeric(),
                                       ID = factor())
                     
                     
                     for(i in 1:length(IPSinput)) {
                       
                       IPSworking <- read.csv(IPSinput[i])
                       IPSworking <- IPSworking[,1:6] # Keep only relevant columns
                       IPSworking$Weighted_Metabolites <- ((IPSworking$Hits.sig + 1) ^ 2) + (IPSworking$Hits.total - IPSworking$Hits.sig)
                       IPSworking$Numerator <- IPSworking$Weighted_Metabolites / (IPSworking$Pathway.total * IPSworking$Expected)
                       IPSworking$Squared_P_Value <- IPSworking$P.Fisher. ^ 2
                       
                       IPSworking$IPS_Value <- IPSworking$Numerator / IPSworking$Squared_P_Value + 1
                       
                       # Clean up (Keep only pathway name, IPS, and Map)
                       IPSworking <- IPSworking[,c(1,10)]
                       IPSworking$log2 <- log2(IPSworking$IPS_Value)
                       names(IPSworking) <- c("Pathway", "IPS", "Log2(IPS)")
                       IPSworking$ID <- gsub(".csv","",input$file3[i,]$name)
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
                     }
                     
                     x <- pheatmap(IPS,
                                   cluster_rows = input$IPSclustRow,
                                   cluster_cols = input$IPSclustCol,
                                   angle_col = 45,
                                   fontsize = 12,
                                   scale = "column",
                     )
                     
                   }, height = 600)
                   
                   output$downloadIPS <- downloadHandler(
                     filename = function(){paste("IPS Value Table.csv")}, 
                     content = function(fname){
                       
                       write.csv(IPS(), fname, row.names = FALSE)
                       
                     }
                   )
                   
                   
                   
                   #Pathway Plots Preview
                   
                   output$Pathways<- renderPlot({
                     req(plotdataframes)
                     ggplot(plotdataframes, aes(x = Enrichment, y = log10Pvalue, label = Label)) + 
                       geom_point(aes(color = log10Pvalue)) + 
                       scale_color_gradient(low = "Yellow", high = "Red") + 
                       ggrepel::geom_text_repel() + 
                       theme_classic()+
                       facet_wrap(~sample, ncol = 1)
                   }
                   )
                   
                   
                   output$AnalysisComplete <- 
                     renderUI({ 
                       show_alert(
                         title= "Analysis Complete", 
                         text= "Please hit the Download Files button to download analysis results", 
                         type="success", 
                         btn_labels="Ok",
                         closeOnClickeOutside=TRUE,
                         ShowCloseButton=TRUE,
                       )
                       
                       
                     })
                   
                   remove_modal_spinner()
                   
                   
                 })
}