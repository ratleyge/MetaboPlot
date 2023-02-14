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
             sidebarPanel(
               
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
               
               h4("P-table & Heatmap"),
               numericInput("featureNumber", "Select top n features for heatmap:", 100),
               actionButton("heatmapButton", "Generate Heatmap"),
               
               h4("Limma"),
               actionButton("limmaButton", "Run Limma"),
               
               
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
               
               h3("NMDS Plot"),
               p("Upload your data and select the ", em("Generate NMDS"), " button."),
               
               # Output: Histogram ----
               
               
               plotOutput(outputId = "NMDS"),
               
               p(),
               
               conditionalPanel(
                 "input.outlierRemove.indexOf('TRUE') > -1",
                 htmlOutput("outliers")
               ),
               
               
               h3("P-table & Heatmap"),
               p("Upload your data and select the ", em("Generate Heatmap"), " button."),
               
               tableOutput("ptable"),
               downloadButton('downloadPtable',"Download P-table"),
               
               plotOutput("heatmap"),
               
               
               h3("Limma"),
               tableOutput("toptable"),
               downloadButton('downloadToptable',"Download Top Table"),
               plotOutput("volcanoPlot"),
               
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
               fileInput("file2", "Upload mummichog pathway enrichment:",
                         multiple = TRUE,
                         accept = c("text/csv",
                                    "text/comma-separated-values,text/plain",
                                    ".csv")),
               
             ),
             
             
             mainPanel(
               
               tableOutput("IPS"),
               downloadButton('downloadIPS',"Download IPS"),
               
             )
           )
           
  )
  
)


