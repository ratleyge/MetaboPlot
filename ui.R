#Fresh
# Define UI for app  ----
ui <- navbarPage(
  header = tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "custom.css")
  ),
  theme = shinytheme("yeti"),
  position = c("fixed-top"),
  title = "MetaboPlot",
  
  # App title ----
  # tabPanel(
  #   "About",
  #   class = "about",
  #   img(
  #     src = 'Logo.jpg',
  #     height = "200px",
  #     width = "100px"
  #   ),
  #   h2("Welcome!"),
  #   p(
  #     "I have made some tools for the ETU lab to process their metabolomics data. Hope they help :)"
  #   )
  # ),
  
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
               
               radioButtons(
                 "metAnnotations",
                 label = "Include annotations:",
                 choices = c(
                   "m/z only" = "mzOnly",
                   "Annotations only" = "annotOnly",
                   "Both" = "mzAnnot"
                 ),
                 selected = "mzOnly"
               ),
               
               h4("Quality Control"),
               numericInput(
                 "missingness",
                 "Remove metabolites with % missingness greater than or equal to:",
                 90,
                 min = 0,
                 max = 100
               ),
               checkboxInput("outlierRemove", "Remove outlier samples.", FALSE),
               
               h4("Groups"),
               uiOutput("Groups"),
               
               fluidRow(id = "relevelRow",
                        
                        column(
                          6,
                          selectInput(
                            inputId = "relevelSelector",
                            label = "Select Control",
                            choices = NULL
                          ),
                        ),
                        
                        column(
                          6, actionButton("relevel", "Relevel Groups")
                        ),),
               actionButton("selectGroups", "Select Groups"),
               textOutput("potato"),
               
               
               h4("Analysis"),
               
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
               
               conditionalPanel(
                 "$.inArray('NMDS', input.Plots) > -1",
                 h4("NMDS settings"),
                 checkboxInput(
                   "ANOSIM",
                   "Calculate ANOSIM statistic for NMDS. This can add a few minutes to the run time.",
                   FALSE
                 ),
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
               downloadButton('downloadPtable', "Download P-tables"),
               h2(" "),
               
               
               tabsetPanel(
                 id = "plots",
                 
                 tabPanel("NMDS", generateNmdsUI("nmdsMod")),
                 tabPanel("Heatmap", generateHeatmapUI("heatmapMod")),
                 tabPanel("Limma", generateLimmaUI("limmaMod")),
                 tabPanel("Lasso", generateLassoUI("lassoMod")),
               ),
             )
           )),
  
  
  
  #Works on its own
  #MetaboAnalyst Tab
  tabPanel("MetaboAnalyst",
           
           sidebarLayout(
             # Sidebar panel for inputs ----
             sidebarPanel(
               id = "MetaboAnalyst form",
               
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
                 "file3",
                 "Upload p-value table as CSV file:",
                 multiple = TRUE,
                 accept = c("text/csv",
                            "text/comma-separated-values,text/plain",
                            ".csv")
               ),
               
               h4("Functional Analysis (Mummichog)"),
               wellPanel(
                 h4("Organisms and Database Options"),
                 radioButtons(
                   "Organisms",
                   label = "Organisms",
                   choices = c(
                     "Homo sapiens (MFN)" = "hsa_mfn",
                     "Homo sapiens (BioCyc)" = "hsa_biocyc",
                     "Homo sapiens (KEGG)" = "hsa_kegg",
                     "Mus musculus (BioCyc)" = "mmu_biocyc",
                     "Mus musculus (KEGG)" = "mmu_kegg",
                     "Escherichia coli K-12 MG1655 (KEGG)" = "eco_kegg"
                   ),
                   
                   selected = "hsa_mfn"
                 ),
                 
                 h4("Ion Mode Options"),
                 radioButtons(
                   "IonMode",
                   label = "Ion Mode",
                   choices = c(
                     "Positive" = "positive",
                     "Negative" = "negative",
                     "Mixed Mode" = "mixed"
                   ),
                   selected = "positive"
                 ),
                 
                 h4("Mass Tolerance (ppm)"),
                 numericInput(
                   "MT",
                   "The instrument's mass accuracy is",
                   5,
                   min = 0,
                   max = 100
                 ),
                 
               ),
               
               actionButton("SubmitMA", "Submit Data"),
               
             ),
             
             # Main panel for displaying outputs ----
             mainPanel(
               h3("Data Preview"),
               p("Upload your data in the panel to the left to preview."),
               
               h3("The Index of Pathway Significance"),
               h4("IPS Table Preview"),
               tableOutput("IPS"),
               downloadButton('downloadIPS', "Download IPS Table"),
               plotOutput("IPSHeatmap"),
               
               h3("Functional Analysis Results"),
               downloadButton('downloadfiles', "Download Files"),
             )
           )),
  
)
