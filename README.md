# MetaboPlot

<p><b>MetaboPlot</b> is an app created for the Epithelial Therapeutics Unit to create simple visualizations and process new data. Currently the app will generate: p-tables, heatmaps, NMDS plots, volcano plots, and lasso feature selection plots. Additionally, we have a section to calculate and visualize the Index of Pathway Significance (IPS).</p>

<h2>Installation and Dependencies</h2>
<p>To run MetaboPlot locally, you'll need to install several dependencies: </p>

``` R
install.packages(pkgs = c("factoextra", "tidyr", "pheatmap", "grid", "ggpubr", "ggrepel", "vegan", "ggrepel", "shinythemes", "glmnet", "glmnetUtilscaret", "BiocManager", "zip", "shinyccsloaders", "caret"), repos='http://cran.rstudio.com/')
# if you are getting an error that limma and Biobase are not available for your version of R, run the line of code below then try again 
# options(BioC_mirror = "http://bioconductor.org")
BiocManager::install("limma")
BiocManager::install("Biobase")
```