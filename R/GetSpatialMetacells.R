#' Construct Spatial Metacells from Single-Cell Spatial Data
#'
#' This function computes spatial metacells from single-cell spatial transcriptomics.
#' Spatial metacells represent aggregated expression profiles for spatially proximal cells.
#' The function utilizes k-nearest neighbor (KNN) weighting to aggregate expression profiles from neighboring cells
#' and generate metacell expression profiles for each cell type.
#'
#' @param normdata A matrix or data frame representing normalized gene expression data,
#' where rows correspond to genes and columns correspond to cells.
#' @param metadata A data frame containing metadata associated with each cell.
#' Must include spatial coordinates (e.g., X and Y) as well as cell type
#' annotations. The row names of the \code{metadata} must match the column names of the \code{normdata}.
#' @param X Character string specifying the column name in \code{metadata} containing the X spatial coordinates.
#' @param Y Character string specifying the column name in \code{metadata} containing the Y spatial coordinates.
#' @param CellType Character string specifying the column name in \code{metadata} containing the cell type annotations.
#' @param spotCoord A data frame containing the spatial coordinates (\code{X} and \code{Y}) of each spatial neighborhood,
#' with neighborhood identifiers as row names.
#' @param k Integer specifying the number of nearest spatial neighbors for constructing metacells.
#' @param radius Numeric value specifying the radius (in units of spatial coordinates) within which neighboring cells are considered.
#' @param min.cells.per.region Integer specifying the minimum number of cells required in each spatial neighborhood to compute metacells.
#' @param ncores Integer specifying the number of CPU cores to use for parallel processing (default: 4).
#'
#' @return A matrix containing the spatial metacell expression profiles, with rows representing individual genes
#'         and columns representing metacells.
#'
#' @examples
#' library(data.table)
#' library(Seurat)
#' library(SpatialEcoTyper)
#' scdata <- fread("https://spatialecotyper.stanford.edu/inc/inc.public.vignettes.php?file=Melanoma1_subset_counts.tsv.gz",
#'                 sep = "\t",header = TRUE, data.table = FALSE)
#' rownames(scdata) <- scdata[, 1]
#' scdata <- as.matrix(scdata[, -1])
#' normdata = NormalizeData(scdata)
#' scmeta <- read.table("https://spatialecotyper.stanford.edu/inc/inc.public.vignettes.php?file=Melanoma1_subset_scmeta.tsv",
#'                       sep = "\t",header = TRUE, row.names = 1)
#'
#' # Construct spatial metacells around B cells
#' snmeta = scmeta[scmeta$CellType=="B", c("X", "Y")]
#' metacells <- GetSpatialMetacells(normdata = normdata,
#'                                  metadata = scmeta,
#'                                  X = "X", Y = "Y",
#'                                  CellType = "CellType",
#'                                  spotCoord = snmeta)
#' head(metacells)
#'
#' @importFrom parallel mclapply
#' @export
#'
GetSpatialMetacells <- function(normdata,
                                metadata,
                                X = "X", Y = "Y",
                                CellType = "CellType",
                                spotCoord = NULL,
                                k = 20, radius = 50,
                                min.cells.per.region = 1,
                                ncores = 4){

  if(!all(c(X, Y, CellType) %in% colnames(metadata))){
    missing_cols = setdiff(c(X, Y, CellType), colnames(metadata))
    stop("Required metadata columns missing: ", paste0(missing_cols, collapse = ", "),
         ". Metadata must include X, Y, and CellType.")
  }
  metadata <- as.data.frame(metadata)
  metadata$CellType <- metadata[, CellType]
  metadata$X <- metadata[, X]
  metadata$Y <- metadata[, Y]
  metadata <- metadata[match(colnames(normdata), rownames(metadata)), ]

  if(!(is.matrix(normdata) | is(normdata, "sparseMatrix"))){
    normdata <- as.matrix(normdata)
  }

  if(is.null(spotCoord)){
    binsize = round(radius*1.4)
    spotCoord = metadata
    spotCoord$SpotID <- paste0("X", round(spotCoord[, X] / binsize),
                              "_Y", round(spotCoord[, Y] / binsize))
    spotCoord <- spotCoord %>% group_by(SpotID) %>%
      summarize(X = median(X), Y = median(Y)) %>% as.data.frame
    rownames(spotCoord) = spotCoord$SpotID
  }else{
    spotCoord$X <- spotCoord[, X]
    spotCoord$Y <- spotCoord[, Y]
  }

  celltypes = table(metadata$CellType)
  celltypes = names(celltypes)[celltypes>k]
  metacell_list <- mclapply(celltypes, function(ct){
    tmpmeta <- metadata[metadata$CellType==ct, ]
    tmpgcm <- normdata[, metadata$CellType==ct]
    weights <- GetKnnWeights(scmeta = tmpmeta, spotCoord,
                             k = k, radius = radius,
                             min.cells.per.region = min.cells.per.region)
    if(is.null(weights)) return(NULL)
    metacell <- tmpgcm %*% weights
    rownames(metacell) <- rownames(normdata)
    colnames(metacell) <- paste0(colnames(weights), "..", ct)
    return(metacell)
  }, mc.cores = ncores)
  names(metacell_list) = celltypes
  metacell_list = metacell_list[lengths(metacell_list)>0]
  ncols = unlist(lapply(metacell_list, ncol))
  metacell_list = metacell_list[ncols>2]
  exclude = setdiff(unique(metadata$CellType), names(metacell_list))
  if(length(exclude)>0) message("Excluding cell types due to insufficient spatial metacells: ",
                                paste0(exclude, collapse = ", "), ".")
  metacell <- do.call(cbind, metacell_list)
  return(metacell)
}

GetKnnWeights <- function(scmeta, spotmeta,
                          k = 20, radius = 50,
                          X = "X", Y = "Y",
                          min.cells.per.region = 1){
  require("Matrix")
  require("RANN")
  require("dplyr")
  locs <- as.data.frame(scmeta[, c(X, Y)])
  colnames(locs) <- c("X", "Y")
  spot_locs = as.data.frame(spotmeta[, c(X, Y)])
  colnames(spot_locs) <- c("X", "Y")
  if(nrow(spot_locs)<3) return(NULL)

  ## sparse distance matrix
  if(k>nrow(locs)) k <- nrow(locs)-1
  knn <- RANN::nn2(data = as.matrix(locs), query = as.matrix(spot_locs), k = k)
  weights <- as(matrix(0, nrow = nrow(locs), ncol = nrow(spot_locs)), "sparseMatrix")
  rownames(weights) <- rownames(locs)
  colnames(weights) <- rownames(spot_locs)

  if(length(weights)>1e8){
    ## Too big integers will be converted to NA
    for(i in 1:ncol(weights)){
      idx = knn$nn.idx[i, ]
      tmp = knn$nn.dists[i, ]+1
      idx2 = idx[idx > 0]
      tmp = tmp[idx > 0]
      idx_tmp = !(is.na(tmp)|is.infinite(tmp))
      tmp = tmp[idx_tmp]
      idx2 = idx2[idx_tmp]
      weights[idx2, i] = tmp
    }
    idx <- Matrix::which(weights>(radius+1), arr.ind = TRUE)
    for(j in unique(idx[,2])){
      i <- idx[idx[,2]==j, 1]
      weights[i, j] <- 0
    }
  }else{
    for(i in 1:k){
      idx <- nrow(weights) * (1:ncol(weights)-1) + knn$nn.idx[, i]
      tmp <- knn$nn.dists[, i]+1
      idx2 <- idx[knn$nn.idx[, i]>0]
      tmp = tmp[knn$nn.idx[, i]>0]
      idx_tmp = !(is.na(tmp)|is.infinite(tmp))
      tmp = tmp[idx_tmp]
      idx2 = idx2[idx_tmp]
      weights[idx2] <- tmp
    }
    weights@x[weights@x>(radius+1)] <- 0
  }
  weights <- drop0(weights)
  ## Unweight
  weights@x <- rep(1, length(weights@x))
  ## Remove spots with fewer than two cells
  weights <- weights[, Matrix::colSums(weights)>=min.cells.per.region]
  if(ncol(weights)<5) return(NULL)
  ## Normalize weights
  weights <- Matrix::t(Matrix::t(weights) / Matrix::colSums(weights))
  return(weights)
}
