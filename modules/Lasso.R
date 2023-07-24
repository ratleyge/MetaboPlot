generateLassoUI <- function (id) {
  
  ns <- NS(id)
  
  tagList(
    h3("Lasso Feature Selection"),
    uiOutput(ns("LassoUIOutput")),
  )
  
}

generateLassoServer <-
  function(id, transdf, plotTitle, groupIdentities) {
    moduleServer(id, function(input, output, session) {
      
      calculateCoefs <- function (transdf) {
        y = transdf$Group
        x = as.matrix(transdf[, 1:(length(transdf)-1)])
        
        # log transform
        x[x == 0] <-1 #Impute zeros
        x <- log2(x) #Log transform
        
        #Pareto Scale 
        x <- paretoScale(x) #apply the function
        
        cv_model <- cv.glmnet(x, 
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
        
        best_model <- glmnet(x, 
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
        return(coefs)
        
      }
      
      # Generate Volcano ----
      
      output$LassoUIOutput <- renderUI({
        FileNumber <- length(levels(transdf$Group))
        plot_output_list_Lasso <- list()
        
        for (i in 2:FileNumber) {
          ns <- session$ns
          
          plot_output_list_Lasso <- append(plot_output_list_Lasso, list(htmlOutput(ns(paste0("LassoTitle", i)))))
          plot_output_list_Lasso <- append(plot_output_list_Lasso, list(tableOutput(ns(paste0("coefs", i)))))
          plot_output_list_Lasso <- append(plot_output_list_Lasso, list(downloadButton(ns(paste0('downloadCoefs', i)), "Download Table")))
          plot_output_list_Lasso <- append(plot_output_list_Lasso, list(
            withSpinner(
              plotOutput(ns(paste0(
                "eachLassoPlot", i
              )), width = "100%", height = "400px"),
              type = getOption("spinner.type", default = 5)
            )
          ))
          
        }
        
        do.call(tagList, plot_output_list_Lasso)
        
      })
      
      
      
      observe({
        FileNumber <- length(levels(transdf$Group))
        
        for (i in 2:FileNumber) {
          local({
            my_i <- i
            tableName <- paste0("coefs", my_i)
            lassoTitle <- paste0("LassoTitle", my_i)
            downloadCoefName <- paste0('downloadCoefs', my_i)
            plotname <- paste0("eachLassoPlot", my_i)
            
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
            workingCoefs <- calculateCoefs(workingTransDf)
            
            output[[lassoTitle]] <-
              renderText(paste0("<h4>", levels(groupIdentities$Group)[my_i], "</h4>"))
            output[[tableName]] <- renderTable(head(workingCoefs))
            
            output[[downloadCoefName]] <- downloadHandler(
              filename = function() {
                paste(plotTitle,
                      "-",
                      levels(groupIdentities$Group)[my_i],
                      "- lasso coefs.csv")
              },
              content = function(fname) {
                write.csv(workingCoefs, fname, row.names = FALSE)
              }
            )
            
            output[[plotname]] <- renderPlot({
              
              workingCoefs <- workingCoefs[abs(workingCoefs$Beta) > 0, ]
              workingCoefs <- workingCoefs[order(-abs(workingCoefs$Beta)), ]
              
              if (class(workingCoefs$m.z) == "numeric") {
                
                g <- (ggplot(workingCoefs, aes(x = reorder(round(m.z, 2), Beta), y = Beta)) +
                        geom_bar(stat = "identity", 
                                 position = position_stack(),
                                 color = "white", 
                                 fill = "lightblue") +
                        ggtitle(plotTitle) +
                        coord_flip() +
                        theme_minimal() + 
                        ylab("Beta Coefficients") +
                        xlab("Input Variables") +
                        theme(plot.margin = unit(c(0.5,1.5,0.5,0.5), "cm")) +
                        theme(text = element_text(size = 15))) +
                  geom_hline(yintercept = 0)
                
              } else {
                
                g <- (ggplot(workingCoefs, aes(x = reorder(m.z, Beta), y = Beta)) +
                        geom_bar(stat = "identity", 
                                 position = position_stack(),
                                 color = "white", 
                                 fill = "lightblue") +
                        ggtitle(plotTitle) +
                        coord_flip() +
                        theme_minimal() + 
                        ylab("Beta Coefficients") +
                        xlab("Input Variables") +
                        theme(plot.margin = unit(c(0.5,1.5,0.5,0.5), "cm")) +
                        theme(text = element_text(size = 15))) +
                  geom_hline(yintercept = 0)
                
              }
              
              g
              
            }, height = 400)
            
            
          })
        }
      })
    })
  }


