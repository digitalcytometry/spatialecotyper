#' Preprocess Spatial Transcriptomics Data
#'
#' This function preprocesses single-cell spatial transcriptomics data by filtering
#' out low-quality genes and cells based on specified thresholds. It ensures that
#' only genes expressed in a minimum number of cells and cells expressing a minimum
#' number of features are retained. Additionally, it reformats the metadata to include
#' spatial coordinates (X and Y).
#'
#' @param expdat A matrix or data frame representing the gene expression data,
#' where rows correspond to genes and columns correspond to cells.
#' @param metadata A data frame containing metadata associated with each cell.
#' Must include spatial coordinates (e.g., X and Y) as well as other cell-specific
#' annotations. The row names of the `metadata` must match the column names of the `expdat`.
#' @param min.cells An integer specifying the minimum number of cells in which a
#' gene must be expressed to be retained (default is 3).
#' @param min.features An integer specifying the minimum number of features (genes)
#' a cell must express to be retained (default is 5).
#' @param X A string specifying the column name in the metadata data frame that
#' represents the X spatial coordinate (default is "X").
#' @param Y A string specifying the column name in the metadata data frame that
#' represents the Y spatial coordinate (default is "Y").
#'
#' @return A list containing two elements:
#' \describe{
#'   \item{expdat}{A filtered matrix of gene expression data, converted to a sparse matrix}
#'   \item{metadata}{A filtered data frame of metadata, aligned with the filtered gene expression data, including reformatted spatial coordinates.}
#' }
#'
#' @examples
#' library(SpatialEcoTyper)
#' library(data.table)
#' library(googledrive)
#' drive_deauth() # no Google sign-in is required
#' drive_download(as_id("1CgUOQKrWY_TG61o5aw7J9LZzE20D6NuI"),
#'                     "HumanMelanomaPatient1_subset_scmeta.tsv", overwrite = TRUE)
#' drive_download(as_id("1CoQmU3u8MoVC8RbLUvTDQmOuJJ703HHB"),
#'               "HumanMelanomaPatient1_subset_counts.tsv.gz", overwrite = TRUE)
#' scdata <- fread("HumanMelanomaPatient1_subset_counts.tsv.gz",
#'                 sep = "\t",header = TRUE, data.table = FALSE)
#' rownames(scdata) <- scdata[, 1]
#' scdata <- as.matrix(scdata[, -1])
#' scmeta <- read.table("HumanMelanomaPatient1_subset_scmeta.tsv",
#'                      sep = "\t", header = TRUE, row.names = 1)
#' processed <- PreprocessST(expdat = scdata, scmeta, X = "X", Y = "Y",
#'                           min.cells = 3, min.features = 5)
#' head(processed$metadata)
#' head(processed$expdat)
#'
#' @export
#'
PreprocessST <- function(expdat, metadata,
                         min.cells = 3,
                         min.features = 5,
                         X = "X", Y = "Y"){
  metadata <- as.data.frame(metadata)
  metadata <- metadata[match(colnames(expdat), rownames(metadata)), ]
  tmp <- sum(rowSums(expdat>0)<min.cells)
  if(tmp>0) message(Sys.time(), " Remove ", tmp, " genes expressed in fewer than ", min.cells, " cells")
  expdat <- expdat[rowSums(expdat>0)>=min.cells, ]
  tmp <- sum(colSums(expdat>0)<min.features)
  if(tmp>0) message(Sys.time(), " Remove ", tmp, " cells with fewer than ", min.features, " features")
  expdat <- expdat[, colSums(expdat>0)>=min.features]

  if(ncol(expdat)>5000) expdat <- as(expdat, "sparseMatrix")
  metadata$X <- metadata[, X]
  metadata$Y <- metadata[, Y]
  metadata <- metadata[match(colnames(expdat), rownames(metadata)), ]
  return(list(expdat = expdat, metadata = metadata))
}

