mystyle <- '
  body {
    padding-top: 70px !important;
  }
'



# Define UI for app  ----
ui <- navbarPage(
  
  header = tags$head(tags$style(type="text/css", mystyle)),
  theme = shinytheme("yeti"),
  position = c("fixed-top"),
  title = "MetaboPlot",
  
  # App title ----
  tabPanel("About", 
      img(src='Logo.jpg', height="200px", width="100px"),
      h2("Welcome!"),
      p("I have made some tools for the ETU lab to process their metabolomics data. Hope they help :)")
    ),
    
  # Sidebar layout with input and output definitions ----
  tabPanel("Analyze",
           
           sidebarLayout(
             
             # Sidebar panel for inputs ----
             sidebarPanel(
               
               # Input: text ----
               textInput("plotTitles", "Type a title for your outputs:", value = "", width = NULL, placeholder = NULL),
               
               # Input: Select a file ----
               fileInput("file1", "Upload feature table:",
                         multiple = TRUE,
                         accept = c("text/csv",
                                    "text/comma-separated-values,text/plain",
                                    ".csv")),
               
               h4("NMDS"),
               checkboxInput("ANOSIM", "Calculate ANOSIM statistic. This can add up to 10 min to the run time.", FALSE),
               actionButton("NMDS", "Generate NMDS plot"),
               
               h4("Heatmaps and P-table"),
               numericInput("featureNumber", "Select top n features for heatmap", 100),
               actionButton("heatmap", "Generate heatmap"),
               
               
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
               p("Upload your data and select the ", em("Generate NMDS plot"), " button."),
               
               # Output: Histogram ----
               plotOutput(outputId = "NMDS"),
               
               tableOutput("outliers")
               
             )
           )
  ), 
  
)


