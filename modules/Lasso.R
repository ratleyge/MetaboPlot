generateLassoUI <- function (id) {
  
  ns <- NS(id)
  
  tagList(
    h3("Lasso Feature Selection"),
    tableOutput(ns("LassoTable")),
    downloadButton(ns('downloadLasso'), "Download Lasso Table"),
    plotOutput(ns("LassoCoefficients")),
  )
  
}


generateLassoServer <- function(id, transdf, plotTitle) {
  
  moduleServer(id, function(input, output, session) {
    
    
    coefs <- isolate({
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
        
      })
    
    
    output$LassoTable <- renderTable({
        
        req(coefs)
        head(coefs)

    }, digits = -1)
    
    
    
    output$downloadLasso <- downloadHandler(
      filename = function(){paste(plotTitles, "- Lasso Coefficients.csv")}, 
      content = function(fname){
        req(coefs)
        write.csv(coefs, fname, row.names = FALSE)
        
      }
    )
    
    output$LassoCoefficients <- renderPlot({
        
        req(coefs)
        
        coefs <- coefs[abs(coefs$Beta) > 0, ]
        coefs <- coefs[order(-abs(coefs$Beta)), ]
        
        g <- (ggplot(coefs, aes(x = reorder(m.z, Beta), y = Beta)) +
                geom_bar(stat = "identity", 
                         position = position_stack(),
                         color = "white", 
                         fill = "lightblue") +
                ggtitle(plotTitle) +
                coord_flip() +
                ylab("Beta Coefficients") +
                xlab("Input Variables") +
                theme(plot.margin = unit(c(0.5,1.5,0.5,0.5), "cm")) +
                theme(text = element_text(size = 15))) +
          geom_hline(yintercept = 0)
        
        g
      
    }, height = length(coefs[abs(coefs$Beta) > 0, ]$Beta) * 20)
    
  })
  
}