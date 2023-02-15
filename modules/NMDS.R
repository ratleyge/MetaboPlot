generateNmdsUI <- function (id) {
  
  ns <- NS(id)
  
  tagList(
    h3("NMDS Plot"),
    plotOutput(outputId = ns("NMDS")),
  )
  
}


generateNmdsServer <- function(id, data.scores, groupIdentities, transdf, anoVal, plotTitle) {
  
  moduleServer(id, function(input, output, session) {
    
    output$NMDS <- renderPlot({
        
        # Add group identities 
        data.scores$Group <- groupIdentities[, 1]
        head(data.scores)
        
        
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
          labs(x = "NMDS1", colour = "Group", y = "NMDS2", shape = "Group", title = paste(plotTitle, "NMDS")) 
        
        
        if (anoVal == TRUE) {
          
          showModal(modalDialog("Calculating ANOSIM...", footer=NULL))
          
          ano <- anosim(transdf, groupIdentities[, 1], distance = "bray", permutations = 9999)
          xx <- xx + labs(subtitle = paste("ANOSIM stat:", 
                                           signif(ano$statistic, digits = 2),
                                           "Significance:", 
                                           signif(ano$signif, digits = 2)))
          
          removeModal()
          
        }
        
        print(xx)

    })
    
  })
  
  
}