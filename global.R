library(factoextra)
library(tidyr)
library(pheatmap)
library(grid)
library(ggpubr)
library(ggrepel)
library(vegan)
library(limma)
library(Biobase)
library(shinythemes)
library(glmnet)
library(glmnetUtils)
library(caret)
library(dplyr)
library(zip)
library(shinycssloaders)
if (require(MetaboAnalystR)) library(MetaboAnalystR)
library(filesstrings)
library(shinybusy)
library(rjson)
library(shinyWidgets)
library(png)
library(gridExtra)
library(fitdistrplus)
library(xlsx)

options(download.file.method = "wininet")

source("modules/NMDS.R")
source("modules/Heatmap.R")
source("modules/Limma.R")
source("modules/Lasso.R")
source("modules/Anova.R")
source("modules/QCfunctions.R")

tempdirect <- paste0(tempdir(), "FunctionalAnalysis", "/")

# Use MIME types only for this list - fileName$type will output only MIME types
# If you add a new type, the logic in server.R also needs to be updated to include it.
acceptedFileTypes <- c("text/csv")
# We will add the following soon
# .xls  - application/vnd.ms-excel
# .xlsx - application/vnd.openxmlformats-officedocument.spreadsheetml.sheet
# .tsv  - text/tab-separated-values

onStop(function() {
  cat("This will run on app stop\n")
  print(ls(envir = .GlobalEnv))
  rm(list = c("transdf"), envir = .GlobalEnv)
  rm(list = c("workingDf"), envir = .GlobalEnv)
  rm(list = c("key"), envir = .GlobalEnv)
  print(ls(envir = .GlobalEnv))
})

    