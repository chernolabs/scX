<h1 align="center">  scXplorer package </h1>
scXplorer is an R package that enables interactive visualization and analysis of single-cell experiment data by creating a Shiny-based application.

## Install
```R
devtools::install_github("tvegawaichman/scXplorer")
```

## Loading Datasets

```R
setwd("/working/directory")
library(devtools)
library(scRNAseq)
load_all()

# Creating SCE object ----
sce <- MarquesBrainData()
cseo <- createSCEobject(xx = sce, 
                        toFactors = "inferred cell type", 
                        toKeep = c("source_name", "age", "Sex", "strain", "treatment"))

launch_scXplorer(cseo)
```
<details><summary> <h2>  Summary </h2> </summary>
  
scXplorer displays a summary of the main descriptive information of the dataset: number of cells and genes, mean number of genes detected per cell, average library size, etc.

![01](/images/01_image_intro.png)

In the summary section, you can explore the relationship between the number of features and the count numbers through graphical visualization.

<p float="left">
  <img src="/images/02_image_summary.png" width="48%" />
  <img src="/images/03_image_summary.png" width="48%" /> 
</p>

</details>

<details><summary> <h2>  Markers </h2> </summary><blockquote>
<details><summary> <h3>  Clusters markers </h3>  </summary><blockquote>
</p>
</blockquote></details>

<details><summary> <h3>  Find new markers </h3> </summary><blockquote>
</p>
</blockquote></details>
</blockquote></details>

## Gene Expression
### Categories
### Fields
### Co-expression
## Differential expression
## Exploratory Data Analysis
### Categories
### Fields
## Visual Tools
### Violin by Partition
### Multiplots

## FAQs

0 - cómo está estructurado un objeto SCE y qué cosas usa el shiny

1 - cómo levantar un objeto 10x, convertirlo en SCE y meterlo en la función

2 - cómo calcular particiones (tipo louvain) y agregarlas al SCE

3 - agregar funcion aparte para el subsampleo (que entre cseo y salga cseo).. y que pongamos que recomendamos hacerlo en el caso de > tantas (100k) células

4 - si tengo un objeto Seurat con un assay específico distinto a los típicos, cómo convertir el objeto Seurat a SCE por fuera de la fución
