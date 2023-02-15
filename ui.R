# Define UI for app  ----
ui <- navbarPage(
  header = tags$head(tags$link(rel = "stylesheet", type = "text/css", href = "custom.css")),
  theme = shinytheme("yeti"),
  position = c("fixed-top"),
  title = "MetaboPlot",
  
  # App title ----
  tabPanel(
    "About",
    class = "about",
    img(
      src = 'Logo.jpg',
      height = "200px",
      width = "100px"
    ),
    h2("Welcome!"),
    p(
      "I have made some tools for the ETU lab to process their metabolomics data. Hope they help :)"
    )
  ),
  
  # Sidebar layout with input and output definitions ----
  tabPanel("Analyze",
           
           sidebarLayout(
             # Sidebar panel for inputs ----
             sidebarPanel(
               id = "form",
               
               h4("Data Upload"),
               
               
               # Input: text ----
               textInput(
                 "plotTitles",
                 "Type a title for your outputs:",
                 value = "",
                 width = NULL,
                 placeholder = NULL
               ),
               
               # Input: Select a file ----
               fileInput(
                 "file1",
                 "Upload feature table:",
                 multiple = TRUE,
                 accept = c("text/csv",
                            "text/comma-separated-values,text/plain",
                            ".csv")
               ),
               checkboxInput("outlierRemove", "Remove outlier samples.", TRUE),
               uiOutput("Groups"),
               selectInput(inputId = "relevelSelector", label = "Select Control", choices = NULL),
               actionButton("relevel", "Relevel groups"),
               p(" "),
               checkboxGroupInput(
                 "Plots",
                 "Analyses to include",
                 c(
                   "NMDS" = "NMDS",
                   "Heatmap" = "heatmap",
                   "Limma & Volcano Plot" = "limma",
                   "Lasso feature selection" = "lasso"
                 ),
               ),
               
               h4("NMDS settings"),
               checkboxInput(
                 "ANOSIM",
                 "Calculate ANOSIM statistic for NMDS. This can add a few minutes to the run time.",
                 FALSE
               ),
               actionButton("Submit", "Submit Data"),
               
             ),
             
             # Main panel for displaying outputs ----
             mainPanel(
               h3("Data Preview"),
               p("Upload your data in the panel to the left to preview."),
               
               tableOutput("viewInputTable"),
               
               conditionalPanel(condition = "input.outlierRemove.indexOf('TRUE') > -1",
                                htmlOutput("outliers")),
               
               h3("P-table"),
               tableOutput("ptable"),
               downloadButton('downloadPtable', "Download P-table"),
               h2(" "),
               
               
               tabsetPanel(id ="plots",
                 
                 tabPanel("NMDS", generateNmdsUI("nmdsMod")),
                 tabPanel("Heatmap", generateHeatmapUI("heatmapMod")),
                 tabPanel("Limma", generateLimmaUI("limmaMod")),
                 tabPanel("Lasso", generateLassoUI("lassoMod")),
               ),
             )
           )),
  
  tabPanel("IPS",
           
           sidebarLayout(
             # Sidebar panel for inputs ----
             sidebarPanel(
               h4("Data Upload"),
               
               # Input: text ----
               #textInput("plotTitles", "Type a title for your outputs:", value = "", width = NULL, placeholder = NULL),
               
               # Input: Select a file ----
               fileInput(
                 "files2",
                 "Upload mummichog pathway enrichment:",
                 multiple = TRUE,
                 accept = c("text/csv",
                            "text/comma-separated-values,text/plain",
                            ".csv")
               ),
               
             ),
             
             
             mainPanel(
               h3("The Index of Pathway Significance"),
               p(
                 paste(
                   "The Index of Pathway Significance (IPS) is an in-house formula
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
                       visualized in figures such as heat maps and metabolic pathway maps."
                 )
               ),
               
               h4("IPS Table Preview"),
               tableOutput("IPS"),
               downloadButton('downloadIPS', "Download IPS Table"),
               
               plotOutput("IPSHeatmap"),
               
             )
           ))
  
)
