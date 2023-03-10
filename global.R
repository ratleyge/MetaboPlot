library(factoextra)
library(tidyr)
library(pheatmap)
library(grid)
library(ggpubr)
library(ggrepel)
library(vegan)
library(limma)
library(Biobase)
library(ggrepel)
library(shinythemes)
library(glmnet)
library(glmnetUtils)
library(caret)

source("modules/NMDS.R")
source("modules/Heatmap.R")
source("modules/Limma.R")
source("modules/Lasso.R")


onStop(function() {
  cat("This will run on app stop\n")
  print(ls(envir = .GlobalEnv))
  rm(list = c("transdf"), envir = .GlobalEnv)
  rm(list = c("workingDf"), envir = .GlobalEnv)
  rm(list = c("key"), envir = .GlobalEnv)
  print(ls(envir = .GlobalEnv))
})

    