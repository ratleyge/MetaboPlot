mystyle <- '
  body {
    padding-top: 70px !important;

  }
  
  #NMDS.shiny-plot-output.shiny-bound-output {
    height: auto !important;
  } 
  
  .navbar {
    font-size: 20px !important;
  }
  
  .navbar-brand {
    font-size: 30px !important;
  }
  
  .navbar-default {
    background-color: #00425A !important;
  }
  
  .about {
    position: absolute;
    top: 50%;
    text-align: center;
    line-height: 2;
    left: 50%;
    transform: translate(-50%, -50%);
  }
  
  .btn {
    color: white !important;
    background-color: #00425A !important;
  }
  
  #heatmap.shiny-plot-output.shiny-bound-output {
    height: auto !important;
  }
  
  #volcanoPlot.shiny-plot-output.shiny-bound-output {
    height: auto !important;
  }
  
'


# Define UI for app  ----
ui <- navbarPage(
  
  header = tags$head(tags$style(type="text/css", mystyle)),
  theme = shinytheme("yeti"),
  position = c("fixed-top"),
  title = "MetaboPlot",
  
  # App title ----
  tabPanel("About", class = "about",
           img(src='Logo.jpg', height="200px", width="100px"),
           h2("Welcome!"),
           p("I have made some tools for the ETU lab to process their metabolomics data. Hope they help :)")
  ),
  
  # Sidebar layout with input and output definitions ----
  tabPanel("Analyze",
           
           sidebarLayout(
             
             # Sidebar panel for inputs ----
             sidebarPanel(id = "form",
               
               h4("Data Upload"),
               
               # Input: text ----
               textInput("plotTitles", "Type a title for your outputs:", value = "", width = NULL, placeholder = NULL),
               
               # Input: Select a file ----
               fileInput("file1", "Upload feature table:",
                         multiple = TRUE,
                         accept = c("text/csv",
                                    "text/comma-separated-values,text/plain",
                                    ".csv")),
               
               h4("NMDS"),
               checkboxInput("outlierRemove", "Remove outlier samples.", TRUE),
               checkboxInput("ANOSIM", "Calculate ANOSIM statistic. This can add up to 10 min to the run time.", FALSE),
               actionButton("NMDS", "Generate NMDS"),
               
               h4("Heatmap"),
               numericInput("featureNumber", "Select top n features for heatmap:", 100),
               actionButton("heatmapButton", "Generate Heatmap"),
               
               h4("Limma & Volcano Plot"),
               actionButton("limmaButton", "Run Limma"),
               
               h4("Lasso"),
               actionButton("lassoButton", "Select Features"),
               
               
             ),
             
             # Main panel for displaying outputs ----
             mainPanel(
               
               h3("Data Preview"),
               p("Upload your data in the panel to the left to preview."),
               
               tableOutput("viewInputTable"),
               
               h3("Group Levels"),
               p("At the moment, this app can only handle data with 2 groups."),
               
               uiOutput("Groups"),
               
               p(),
               
               actionButton("relevel", "Relevel groups"),
               
               h3("P-table"),
               tableOutput("ptable"),
               downloadButton('downloadPtable',"Download P-table"),
               
               h3("NMDS Plot"),
               p("Upload your data and select the ", em("Generate NMDS"), " button."),
               
               # Output: Histogram ----
               
               
               plotOutput(outputId = "NMDS"),
               
               p(),
               
               conditionalPanel(
                 "input.outlierRemove.indexOf('TRUE') > -1",
                 htmlOutput("outliers")
               ),
               
               
               h3("Heatmap"),
               p("Upload your data and select the ", em("Generate Heatmap"), " button."),
               plotOutput("heatmap"),
               
               
               h3("Limma & Volcano Plot"),
               p("Upload your data and select the ", em("Run Limma"), " button."),
               tableOutput("toptable"),
               downloadButton('downloadToptable',"Download Top Table"),
               plotOutput("volcanoPlot"),
               
               h3("Lasso Feature Selection"),
               p("Upload your data and select the ", em("Select Features"), " button."),
               tableOutput("LassoTable"),
               downloadButton('downloadLasso',"Download Lasso Table"),
               plotOutput("LassoCoefficients"),

               
             )
           )
  ), 
  
  tabPanel("IPS", 
           
           sidebarLayout(
             
             # Sidebar panel for inputs ----
             sidebarPanel(
               
               h4("Data Upload"),
               
               # Input: text ----
               #textInput("plotTitles", "Type a title for your outputs:", value = "", width = NULL, placeholder = NULL),
               
               # Input: Select a file ----
               fileInput("files2", "Upload mummichog pathway enrichment:",
                         multiple = TRUE,
                         accept = c("text/csv",
                                    "text/comma-separated-values,text/plain",
                                    ".csv")),
               
             ),
             
             
             mainPanel(
               
               h3("The Index of Pathway Significance"),
               p(paste("The Index of Pathway Significance (IPS) is an in-house formula
                       derived from five metrics provided by Metaboanalyst: 
                       total number of metabolites in a metabolic pathway, the total number 
                       of metabolite hits picked up from a sample, the number of significant 
                       metabolite hits in a pathway, the expected number of total metabolite 
                       hits in a sample, and FET, the p-value associated with a sample. 
                       When considered together, they produce a single, easily comparable value 
                       that is associated with differences between two treatments with respect to 
                       a certain metabolic pathway. A higher IPS value indicates a greater difference 
                       and can be used to rank either the metabolic pathways that are shared between two
                       treatments, or the treatments that all exhibit that particular pathway. Assigning 
                       a single value to each pathway also allows many comparisons to be condensed and 
                       visualized in figures such as heat maps and metabolic pathway maps.")),
               
               h4("IPS Table Preview"),
               tableOutput("IPS"),
               downloadButton('downloadIPS',"Download IPS Table"),
               
               plotOutput("IPSHeatmap"),
               
             )
           )
           
  )
  
)


