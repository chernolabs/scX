<h1 align="center">  scXplorer package </h1>

scXplorer is an R package that enables interactive visualization and analysis of single-cell experiment data by creating a Shiny-based application. With scXplorer all aspects of single-cell data can be explored using a large number of different types of plots, such as scatter plots, heatmaps, boxplots, dot and violins plots, etc. All the information associated with the cells can be displayed in a customized way: both numerical variables such as logcounts or pseudotime values, and categorical variables such as cell types or sample. One of the main hallmarks of scXplorer is the possibility to plot the main embeddings used for single cell - UMAP, tSNE and PCA - both in 2D and 3D in an interactive way. Thus, embeddings can be rotated, translated and zoomed. But scXplorer is not only a visualization tool, it also allows you to perform different types of analysis on single cell data, such as finding the markers of a cell type or determining the differential genes between two different conditions.

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
library(scRNAseq)
load_all()
sce <- MarquesBrainData()
```
scXplorer app can be created and launched with only two functions. createSCEobject() creates the single cell experiment object that will be used inside the app. With this function we can configure the main features of our dataset, such as what is the main partition and what metadata information we want to analyze with scXplorer.  On the other hand, different aspects of the data preprocessing are also defined. Finally with launch_scXplorer() we can initiate the app.

```R
# Creating SCE object
cseo <- createSCEobject(xx = sce, 
                        toFactors = "inferred cell type", 
                        toKeep = c("source_name", "age", "Sex", "strain", "treatment"))

launch_scXplorer(cseo)
```
## Summary
  
scXplorer displays a summary of the main descriptive information of the dataset: number of cells and genes, mean number of genes detected per cell, average library size, etc.

![01](/images/01_image_intro.png)

In the summary section, you can explore the relationship between the number of features and the count numbers through graphical visualization.

<p float="left">
  <img src="/images/02_image_summary.png" width="48%" />
  <img src="/images/03_image_summary.png" width="48%" /> 
</p>

</details>

##  Markers

In "Markers" section there are two types of analysis. On the one hand, in "Cluster markers" clicking on a cell displays a table with the marker genes of the cluster to which that cell belongs. On the other hand, in "Find new markers" you can select a set of cells in the embedding and scXplorer will calculate their marker genes.

<details><summary> <h3>  Clusters markers </h3>  </summary><blockquote>
 
 Aca tenes el gif
 
<img src="/images/scXplorer-Google-Chrome-2023-07-05-22-45-05.gif" width="80%" />

</blockquote></details>

<details><summary> <h3>  Find new markers </h3> </summary><blockquote>
</p>
</blockquote></details>
</blockquote></details>

##  Gene Expression
<details><summary> <h3> Categories </h3>  </summary><blockquote>
</p>
</blockquote></details>

<details><summary> <h3>  Fields </h3> </summary><blockquote>
</p>
</blockquote></details>
</blockquote></details>

<details><summary> <h3>  Co-expression </h3> </summary><blockquote>
</p>
</blockquote></details>
</blockquote></details>

## Differential expression

##  Exploratory Data Analysis
<details><summary> <h3> Categories </h3>  </summary><blockquote>
</p>
</blockquote></details>

<details><summary> <h3>  Fields </h3> </summary><blockquote>
</p>
</blockquote></details>
</blockquote></details>

##  Visual Tools
<details><summary> <h3> Violin by Partition </h3>  </summary><blockquote>
</p>
</blockquote></details>

<details><summary> <h3>  Multiplots </h3> </summary><blockquote>
</p>
</blockquote></details>
</blockquote></details>

## FAQs

0 - cómo está estructurado un objeto SCE y qué cosas usa el shiny

1 - cómo levantar un objeto 10x, convertirlo en SCE y meterlo en la función

2 - cómo calcular particiones (tipo louvain) y agregarlas al SCE

3 - agregar funcion aparte para el subsampleo (que entre cseo y salga cseo).. y que pongamos que recomendamos hacerlo en el caso de > tantas (100k) células

4 - si tengo un objeto Seurat con un assay específico distinto a los típicos, cómo convertir el objeto Seurat a SCE por fuera de la fución
