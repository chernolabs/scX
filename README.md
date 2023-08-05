<h1 align="center">  scXplorer package </h1>

scXplorer is an R package that enables interactive visualization and analysis of single-cell experiment data by creating a Shiny-based application. With scXplorer all aspects of single-cell data can be explored using a large number of different types of plots, such as scatter plots, heatmaps, boxplots, dot and violins plots, etc. All the information associated with cells can be displayed in a customized way: both numerical variables such as logcounts or pseudotime values, and categorical variables such as cell types or sample. One of the main hallmarks of scXplorer is the possibility to plot the main embeddings used for single cell - UMAP, tSNE and PCA - both in 2D and 3D in an interactive way. Thus, embeddings can be rotated, translated and zoomed. But scXplorer is not only a visualization tool, it also allows you to perform different types of analysis on single cell data, such as finding the markers of a cell type or determining the differential genes between two different conditions.

## Table of Contents
- [Install](#install)
- [Quick Start Guide](#quick-start-guide)
  * [Summary](#summary)
  * [Markers](#markers)
  * [Gene expression](#gene-expression)
  * [Differential expression](#differential-expression)
  * [Exploratory Data Analysis](#exploratory-data-analysis)
  * [Visual Tools](#visual-tools)
- [FAQs](#faqs)


## Install
scXplorer can be installed from Github as follows:
```R
devtools::install_github("tvegawaichman/scXplorer")
```
## Quick Start Guide
### Loading Datasets
To show the different features that scXplorer has we will use single cell data related to the oligodendrocyte developmental lineage ([Marques et al. 2016](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5221728/)). In this dataset we have cells from 12 regions of the central nervous system of juvenile and adult mice, and 10 distinct cell populations were identified. The data can be obtained using the "scRNAseq" package as follows:

```R
setwd("/working/directory")
library(devtools)
load_all(quiet=TRUE)
sce <- readRDS("./data/ssce.RDS")
```
scXplorer app can be created and launched with only two functions. createSCEobject() creates the single cell experiment object that will be used inside the app. With this function we can configure the main features of our dataset, such as what is the main partition and what metadata information we want to analyze with scXplorer.  On the other hand, different aspects of the data preprocessing are also defined. Finally with launch_scXplorer() we can initiate the app.

```R
# Creating SCE object
cseo <- createSCEobject(xx = sce, 
                        toFactors = "inferred_cell_type", 
                        toKeep = c("source_name", "age", "sex", "strain", "treatment"),
                        descriptionText = "Quick Start Guide")

launch_scXplorer(cseo)
```
## Summary
  
scXplorer displays a summary of the main descriptive information of the dataset: number of cells and genes, mean number of genes detected per cell, average library size, etc.

In the summary section, you can explore the relationship between the number of features and the count numbers through graphical visualization.

<img src="/images/summary.gif" width="100%" />

</details>

##  Markers

In "Markers" section there are two types of analysis. On the one hand, in "Cluster markers" clicking on a cell displays a table with the marker genes of the cluster to which that cell belongs. On the other hand, in "Find new markers" you can select a group of cells in the embedding and scXplorer will calculate their marker genes.

<details><summary> <h3>  Clusters markers </h3>  </summary><blockquote>

This section allows you to find the marker genes for the partition defined in the single-cell object, typically cell types or cell states. Clicking on one of the cells in the embedding will display a table of marker genes for the partition to which that cell belongs. For each of the markers different metrics such as boxcor, robustness and FDR are displayed. This table can be downloaded in various formats, such as .csv, .xlsx .pdf, or you can copy it to the clipboard. 

By clicking on a marker in the table you can see its expression profile across the entire dataset in the embedding. In addition violin and spikeplots are displayed at the bottom.
 
<img src="/images/cluster_markers.gif" width="100%" />

</blockquote></details>

<details><summary> <h3>  Find new markers </h3> </summary><blockquote>

Here you can select with the box or lasso tool a set of cells in the embedding and scXplorer will calculate the marker genes. You can download not only the marker table but also the selected cell list.

As in the previous section, if you click on one of the markers you can see its expression along the dataset with violin and spikeplots.

<img src="/images/new_markers.gif" width="100%" />

</blockquote></details>
</blockquote></details>

##  Gene Expression

In "Gene Expression" you can explore different aspects of the expression of one or more genes of interest. Determine how expression changes according to different categorical and continuous variables, as well as analyse co-detection between pairs of genes.

<details><summary> <h3> Categories </h3>  </summary><blockquote>

In Settings you can select one or more genes or upload a file with a list of genes. The average expression of the genes of interest can be viewed in the different embeddings available, with the possibility to colour according to the different SCE partitions to compare gene expression with different cell types or conditions present in the metadata. 

A wide variety of plots are available to analyse the expression of the genes of interest in the different categories. Heatmaps allow normalisation of expression by gene, clustering by row and column and grouping of cells by condition. Similarly, dotplots allow normalisation of expression and clustering of genes. 

<img src="/images/categories.gif" width="100%" />

</blockquote></details>

<details><summary> <h3>  Fields </h3> </summary><blockquote>

**Fields** allows you to analyse the expression of a set of genes in relation to numeric variables present in your dataset, such as the number of counts or pseudotime value, if present in the metadata of the sce object. Below the embedding, a line plot of the average expression of the genes of interest as a function of the chosen variables and a spikeplot are displayed. Furthermore, you can find heatmaps sorted by the chosen numerical variable that can be divided according to some categorical variable, and multiline plots showing the comparison of the expression profile of the genes of interest along the field.

<img src="/images/ge_fields.gif" width="100%" />

</blockquote></details>
</blockquote></details>

<details><summary> <h3>  Co-expression </h3> </summary><blockquote>

In **Co-expression** section you can analyse the co-appearance of pairs of genes, determine the number and percentage of cells in which each gene is expressed separately and together. You can also view this information graphically in the embedding. In addition, the co-expression of these genes in the different conditions of any of the partitions in the dataset can be analysed by a co-detection matrix.

<img src="/images/coexpression.gif" width="100%" />

</blockquote></details>
</blockquote></details>

## Differential expression

Coming soon ...
<img src="/images/differential_expression.gif.gif" width="100%" />

##  Exploratory Data Analysis
<details><summary> <h3> Categories </h3>  </summary><blockquote>

<img src="/images/exploratory_categories.gif" width="100%" />

</blockquote></details>

<details><summary> <h3>  Fields </h3> </summary><blockquote>

<img src="/images/fields.gif" width="100%" />

</blockquote></details>
</blockquote></details>

##  Visual Tools
<details><summary> <h3> Violin by Partition </h3>  </summary><blockquote>

<img src="/images/violins.gif" width="100%" />

</blockquote></details>

<details><summary> <h3>  Multiplots </h3> </summary><blockquote>

<img src="/images/multiplots.gif" width="100%" />

</blockquote></details>
</blockquote></details>

## FAQs

0 - cómo está estructurado un objeto SCE y qué cosas usa el shiny

1 - cómo levantar un objeto 10x, convertirlo en SCE y meterlo en la función

2 - cómo calcular particiones (tipo louvain) y agregarlas al SCE

3 - agregar funcion aparte para el subsampleo (que entre cseo y salga cseo).. y que pongamos que recomendamos hacerlo en el caso de > tantas (100k) células

4 - si tengo un objeto Seurat con un assay específico distinto a los típicos, cómo convertir el objeto Seurat a SCE por fuera de la fución
