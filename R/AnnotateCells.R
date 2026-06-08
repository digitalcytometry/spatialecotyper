#' Extract Spatial Ecotype Annotations for Single Cells
#'
#' This function adds spatial ecotype annotations to the metadata of single cells.
#'
#' @param scmeta Data frame containing single-cell metadata. Two columns (X and Y) for spatial coordinates are required.
#' @param obj A Seurat object returned from SpatialEcoTyper function.
#' @param col Character string specifying the name of the column in `obj@meta.data`
#' containing the spatial ecotype annotations. Default is `"SE"`.
#' @param dropcell Logical. If TRUE, cells that cannot be assigned to any spatial
#' ecotype (outside the radius) will be removed from the returned metadata. Default is TRUE.
#' @return An updated version of `scmeta` with spatial ecotype annotations added.
#'
#' @examples
#' library(data.table)
#' library(Seurat)
#' library(SpatialEcoTyper)
#'
#' scmeta <- read.table("https://spatialecotyper.stanford.edu/inc/inc.public.vignettes.php?file=Melanoma1_subset_scmeta.tsv",
#'                      sep = "\t", header = TRUE, row.names = 1)
#' head(scmeta)
#'
#' url <- "https://spatialecotyper.stanford.edu/inc/inc.public.vignettes.php?file=Melanoma1_subset_SpatialEcoTyper_results.rds"
#' download.file(url, destfile = "Melanoma1_subset_SpatialEcoTyper_results.rds", mode = "wb")
#' obj <- readRDS("Melanoma1_subset_SpatialEcoTyper_results.rds")$obj
#'
#' ## Transfer SE annotations to single cells
#' scmeta <- AnnotateCells(scmeta = scmeta, obj = obj, dropcell = TRUE)
#' head(scmeta)
#'
#' @export
#'
AnnotateCells <- function(scmeta, obj, col = "SE", dropcell = TRUE){
  if(!all(c("X", "Y") %in% colnames(scmeta))){
    stop("Required columns missing from scmeta: ", paste0(c("X", "Y"), collapse = ", "), ". Please include X and Y coordinates.")
  }
  if(!all(c("X", "Y") %in% colnames(obj@meta.data))){
    stop("The provided Seurat object does not appear to be a valid SpatialEcoTyper result. It must include 'X' and 'Y' in obj@meta.data.")
  }
  radius = as.numeric(gsub("_.*", "", gsub(".*radius", "", obj@project.name)))
  ### One nearest neighbor
  knn <- RANN::nn2(data = as.matrix(obj@meta.data[, c("X", "Y")]),
                   query = as.matrix(scmeta[, c("X", "Y")]), k = 1)
  scmeta[, col] <- obj@meta.data[knn$nn.idx[,1], col]
  scmeta[knn$nn.dists[,1]>radius, col] <- NA
  if(dropcell) scmeta = scmeta[!is.na(scmeta[, col]), ]
  return(scmeta)
}
