#' @title SCE dataset from Marques et al. 2016 
#'
#' @description A subset of data obtained from \href{https://bioconductor.org/packages/release/data/experiment/html/scRNAseq.html}{scRNAseq package}. 
#' In this dataset from \href{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5221728/}{Marques et al. 2016} we have cells from 12 regions of the central nervous system of juvenile and adult mice and 10 distinct cell populations 
#' have been identified. This package includes a modified version of this dataset in which the original cells have been subsampled and a pseudotime has been calculated to demonstrate how scXplorer can represent numerical attributes. 
#' Pseudotime was calculated with \href{https://www.bioconductor.org/packages/release/bioc/html/slingshot.html}{slingshot}. 
#'
#' @format ## `sce`
#' A Single Cell Experiment object with 2521 cells and 8691 genes:
#' \describe{
#'   \item{assays}{counts}
#'   \item{colData}{title source_name age inferred_cell_type sex strain treatment pseudotime}
#' }
#' @source <https://bioconductor.org/packages/release/data/experiment/manuals/scRNAseq/man/scRNAseq.pdf#Rfn.SingleCellExperiment.Rdash.class>
"sce"