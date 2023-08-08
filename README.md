<h1 align="center">  scXplorer package </h1>


<div align="justify">
<img align="left" width="40%" src="inst/www/scXplorer-03.png"> scXplorer is an R package that enables interactive visualization and analysis of single-cell experiment data by creating a Shiny-based application. With scXplorer all aspects of single-cell data can be explored using a large number of different types of plots, such as scatter plots, heatmaps, boxplots, dot and violins plots, etc. All the information associated with cells can be displayed in a customized way: both numerical variables such as logcounts or pseudotime values, and categorical variables such as cell types or sample. One of the main hallmarks of scXplorer is the possibility to plot the main embeddings used for single cell - UMAP, tSNE and PCA - both in 2D and 3D in an interactive way. Thus, embeddings can be rotated, translated and zoomed. But scXplorer is not only a visualization tool, it also allows you to perform different types of analysis on single cell data, such as finding the markers of a cell type or determining the differential genes between two different conditions.
</div>

<h1 align="center">   </h1>

<table>
<tr ><td>
<h3 align="left">  Table of Contents </h3>

- [Install](#install)
- [Quick Start Guide](#quick-start-guide)
  * [Summary](#summary)
  * [Markers](#markers)
  * [Gene expression](#gene-expression)
  * [Differential expression](#differential-expression)
  * <nobr> [Exploratory Data Analysis](#exploratory-data-analysis) </nobr>
  * [Visual Tools](#visual-tools)
- [FAQs](#faqs)
</td>
<td>
<img src="/images/summary.gif" width="100%" />
</td>
</tr>
</table>

## Install
scXplorer can be installed from Github as follows:
```R
devtools::install_github("tvegawaichman/scXplorer")
```
## Quick Start Guide
### Loading Datasets


To show the different features that scXplorer has we will use single cell data related to the oligodendrocyte developmental lineage ([Marques et al. 2016](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5221728/)). In this dataset we have cells from 12 regions of the central nervous system of juvenile and adult mice, and 10 distinct cell populations were identified. This package includes a modified version of this dataset in which the original cells have been subsampled and a pseudotime has been calculated to show how scXplorer can represent numerical attributes.


```R
setwd("/working/directory")
load_all(quiet=TRUE)
sce <- readRDS("scXplorer/data/ssce.RDS")
```

<p align="justify">  
The scXplorer app can be created and launched with only two functions. createSCEobject() creates the single cell experiment object that will be used within the application. With this function we can configure the main features of our dataset, such as which is the main partition and which metadata information we want to analyse with scXplorer.  On the other hand, various aspects of data pre-processing are also defined. Finally we can launch the application with launch_scXplorer().
</p>

```R
# Creating SCE object
cseo <- createSCEobject(xx = sce, 
                        toFactors = "inferred_cell_type", 
                        toKeep = c("source_name", "age", "sex", "strain", "treatment", "pseudotime"),
                        descriptionText = "Quick Start Guide")

launch_scXplorer(cseo)
```
## Summary

<p align="justify">  
scXplorer displays a summary of the main descriptive information of the dataset: number of cells and genes, mean number of genes detected per cell, average library size, etc.

In the summary section, you can explore the relationship between the number of features and the count numbers through graphical visualization.
</p>

<details><summary> ... </summary><blockquote>
<img src="/images/summary.gif" width="100%" />
</blockquote></details>

##  Markers

<p align="justify">  
In "Markers" section there are two types of analysis. On the one hand, in "Cluster markers" clicking on a cell displays a table with the marker genes of the cluster to which that cell belongs. On the other hand, in "Find new markers" you can select a group of cells in the embedding and scXplorer will calculate their marker genes.
</p>

<details><summary> ... </summary><blockquote>
 
###  Clusters markers 

<p align="justify">  
This section allows you to find the marker genes for the partition defined in the single-cell object, typically cell types or cell states. Clicking on one of the cells in the embedding will display a table of marker genes for the partition to which that cell belongs. For each of the markers different metrics such as boxcor, robustness and FDR are displayed. This table can be downloaded in various formats, such as .csv, .xlsx .pdf, or you can copy it to the clipboard. 

By clicking on a marker in the table you can see its expression profile across the entire dataset in the embedding. In addition violin and spikeplots are displayed at the bottom.
</p>


<img src="/images/cluster_markers.gif" width="100%" />


### Find new markers 

<p align="justify">  
Here you can select with the box or lasso tool a set of cells in the embedding and scXplorer will calculate the marker genes. You can download not only the marker table but also the selected cell list.

As in the previous section, if you click on one of the markers you can see its expression along the dataset with violin and spikeplots.
</p>


<img src="/images/new_markers.gif" width="100%" />
</blockquote></details>


##  Gene Expression

<p align="justify">  
In "Gene Expression" you can explore different aspects of the expression of one or more genes of interest. Determine how expression changes according to different categorical and continuous variables, as well as analyse co-detection between pairs of genes.
</p>

<details><summary> ... </summary><blockquote>
 
### Categories

<p align="justify"> 
In Settings you can select one or more genes or upload a file with a list of genes. The average expression of the genes of interest can be viewed in the different embeddings available, with the possibility to colour according to the different SCE partitions to compare gene expression with different cell types or conditions present in the metadata. 

A wide variety of plots are available to analyse the expression of the genes of interest in the different categories. Heatmaps allow normalisation of expression by gene, clustering by row and column and grouping of cells by condition. Similarly, dotplots allow normalisation of expression and clustering of genes. 
</p>


<img src="/images/categories.gif" width="100%" />


### Fields

<p align="justify"> 
Fields allows you to analyse the expression of a set of genes in relation to numeric variables present in your dataset, such as the number of counts or pseudotime value, if present in the metadata of the sce object. Below the embedding, a line plot of the average expression of the genes of interest as a function of the chosen variables and a spikeplot are displayed. Furthermore, you can find heatmaps sorted by the chosen numerical variable that can be divided according to some categorical variable, and multiline plots showing the comparison of the expression profile of the genes of interest along the field.
</p>


<img src="/images/ge_fields.gif" width="100%" />



### Co-expression

<p align="justify"> 
In Co-expression section you can analyse the co-appearance of pairs of genes, determine the number and percentage of cells in which each gene is expressed separately and together. You can also view this information graphically in the embedding. In addition, the co-expression of these genes in the different conditions of any of the partitions in the dataset can be analysed by a co-detection matrix.
</p>


<img src="/images/coexpression.gif" width="100%" />

</blockquote></details>

## Differential expression
<p align="justify"> 
In this section, the differentially expressed genes between two clusters can be obtained. The algorithm used is `findMarkers()` from the scran package and the test type can be defined in the `createSCEobject()` function (see "Summary" section). Different significance levels can be chosen interactively for both the logFC and the FDR. Once the clusters have been chosen, a table will be displayed with the list of differentially expressed genes accompanied by the logFC and the FDR obtained for each of them. This table can be downloaded in different formats: csv, pdf, xlsx. 

The main figure in this section is a VolcanoPlot in which genes that are down and up expressed are coloured in blue and red, respectively. On the other hand, you can also display ViolinPlots, Spikeplots, Heatmaps and Dotplots of the differentially expressed genes.
</p>

<details><summary> ... </summary><blockquote>
<img src="/images/differential_expression.gif" width="100%" />
</blockquote></details>

##  Exploratory Data Analysis

In **Exploratory Data Analysis** section you will be able to understand the relationship between different features contained as metadata in your SCE object.

<details><summary> View more </summary><blockquote>
 
### Categories

<p align="justify"> 
Here you can observe the proportion of cells belonging to the different levels of a categorical variable presented in the metadata and disaggregate these proportions according to the levels of another categorical variable. All this information is displayed in the form of a barplot. In the subsection "Matrix" a confusion matrix between the two selected features can be plotted with the option to display the Jaccard index for each of the grid cells. In addition, the Rand index is displayed, which is a global measure of the similarity between the two clusterings.
</p>

<img src="/images/exploratory_categories.gif" width="100%" />



### Fields

In a similar way to the previous subsection, in "Field" you can explore how the value of one or more numerical variables changes as a function of another variable, either numerical or categorical. You can make different types of plots such as: Distribution Plots, Heatmaps, Dotplots and StackedViolins.

<img src="/images/fields.gif" width="100%" />


</blockquote></details>

##  Visual Tools

In the "Visual Tools" section, different plots can be obtained to explore in more depth different aspects of the single-cell experiment data, such as how the expression of a given set of genes varies at different levels of a feature or to recognise a set of cells of interest within an embedding.

<details><summary> ... </summary><blockquote>

### Violin by Partition 

By selecting a set of genes of interest, a set of ViolinPlots can be obtained for each gene showing its expression at different levels of a feature. These plots can be divided according to the levels of another categorical feature.

<img src="/images/violins.gif" width="100%" />

###  Multiplots

 Multiplots allows you to explore how different variables change across cells in an embedding of your choice, such as the expression of a given set of genes, the partitions of a categorical variable or the value of a continuous variable.
 
<img src="/images/multiplots.gif" width="100%" />


</blockquote></details>

## FAQs

<details><summary> <h3> 1 -   What type of input does scXplorer accept? </h3> </summary><blockquote>

</blockquote></details>

<details><summary> <h3> 2 -  Is there a maximum input size that scXplorer can receive? </h3> </summary><blockquote>

</blockquote></details>

<details><summary> <h3> 3 -  How can I deploy my scXplorer app online? </h3> </summary><blockquote>

</blockquote></details>

