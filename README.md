# MetaboPlot

<p><b>MetaboPlot</b> is an app created for the Epithelial Therapeutics Unit to create simple visualizations and process new data. Currently the app will generate: p-tables, heatmaps, NMDS plots, volcano plots, and lasso feature selection plots. Additionally, we have a section to calculate and visualize the Index of Pathway Significance (IPS).</p>

<h2>Installation and Dependencies</h2>
<p>To run MetaboPlot, you'll need to install several dependencies: </p>

``` R
install.packages(
  pkgs = c("factoextra", "tidyr", "pheatmap", "grid", "ggpubr", "ggrepel", "vegan", "limma", "Biobase",
  "shinythemes","glmnet", "glmnetUtils", "caret", "dplyr", "zip", "shinycssloaders", "MetaboAnalystR",
  "filesstrings", "shinybusy", "rjson", "shinyWidgets", "png", "gridExtra", "fitdistrplus", "BiocManager"), 
  repos = 'http://cran.rstudio.com/'
)
BiocManager::install("limma")
BiocManager::install("Biobase")

# if you are getting an error that limma and Biobase are not available for your version of R, run the line of code below then try again 
options(BioC_mirror = "http://bioconductor.org")
```

<p>If you wish to do pathway analysis using MetaboAnalyst, you'll also need to install MetaboAnalystR. <a html = "https://github.com/xia-lab/MetaboAnalystR">See their repository<a/> for installation instructions: https://github.com/xia-lab/MetaboAnalystR</p>

<h2>Data Format</h2>
<p>Metaboplot accepts metabolomics data in the following format:</p>
<ul>
<li>Columns should be samples</li>
<li>Rows should be metabolites</li>
<li>The first row should contain group identities</li>
<li>The column containing annotations should be labeled "Name"</li>
<li>The column containing m/z values should be labeled "Bucket.label"</li>
</ul>
<p>An example data set can be found in the repository. Coagulase negative staph was exposed to several topical products and MALDI timsTOF mass spectrometry was performed. Samples were annotated with the Collision Cross-Section Consortium database in MetaboScape.</p>










