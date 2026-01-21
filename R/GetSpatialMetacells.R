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
#' @param k Integer specifying the number of nearest spatial neighbors for constructing metacells.
#' @param radius Numeric value specifying the radius (in units of spatial coordinates) within which neighboring cells are considered.
#' @param bin Character string specifying the column name in \code{metadata} containing the spatial neighborhood identifier.
#' @param bin.X Character string specifying the column name in \code{metadata} containing the X spatial coordinates of spatial neighborhoods.
#' @param bin.Y Character string specifying the column name in \code{metadata} containing the Y spatial coordinates of spatial neighborhoods.
#' @param min.cells.per.region Integer specifying the minimum number of cells required in each spatial neighborhood to compute metacells.
#' @param ncores Integer specifying the number of CPU cores to use for parallel processing (default: 1).
#'
#' @return A matrix containing the spatial metacell expression profiles, with rows representing individual genes
#'         and columns representing metacells.
#'
#' @examples
#' library(data.table)
#' library(Seurat)
#' library(SpatialEcoTyper)
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
#' normdata = NormalizeData(scdata)
#' scmeta <- read.table("HumanMelanomaPatient1_subset_scmeta.tsv",
#'                       sep = "\t",header = TRUE, row.names = 1)
#'
#' # Construct spatial metacells from single-cell spatial data
#' metacells <- GetSpatialMetacells(normdata = normdata, metadata = scmeta,
#'                                  X = "X", Y = "Y", CellType = "CellType")
#' head(metacells)
#'
#' @importFrom parallel detectCores mclapply
#' @export
#'

GetSpatialMetacells <- function(normdata, metadata,
                                X = "X", Y = "Y",
                                CellType = "CellType",
                                k = 20, radius = 50,
                                bin = "SpotID",
                                bin.X = "Spot.X",
                                bin.Y = "Spot.Y",
                                min.cells.per.region = 1,
                                ncores = 1){

  metadata <- as.data.frame(metadata)
  metadata <- metadata[match(colnames(normdata), rownames(metadata)), ]

  if(!(is.matrix(normdata) | is(normdata, "sparseMatrix"))){
    normdata <- as.matrix(normdata)
  }
  metadata$CellType <- metadata[, CellType]
  if(!all(c(bin, bin.X, bin.Y) %in% colnames(metadata))){
    metadata$X = metadata[, X]
    metadata$Y = metadata[, Y]
    binsize = round(radius*1.4)
    metadata$SpotID <- paste0("X", round(metadata[, X] / binsize),
                              "_Y", round(metadata[, Y] / binsize))
    metadata <- metadata %>% group_by(SpotID) %>%
      mutate(Spot.X = median(X), Spot.Y = median(Y)) %>% as.data.frame
    bin.X = "Spot.X"
    bin.Y = "Spot.Y"
    bin = "SpotID"
  }
  metadata$SpotID <- metadata[, bin]
  celltypes <- metadata %>% dplyr::count(CellType, SpotID) %>%
    dplyr::count(CellType) %>% filter(n>4) %>% pull(CellType)
  exclude <- setdiff(unique(metadata$CellType), celltypes)
  if(length(exclude)>0) message("\t\tExclude ", paste0(exclude, collapse = ", "), " due to limited number of spatial metacells")

  metacell_list <- mclapply(celltypes, function(ct){
    tmpmeta <- metadata[metadata$CellType==ct, ]
    tmpgcm <- normdata[, metadata$CellType==ct]
    weights <- GetKnnWeights(metadata = tmpmeta, k = k, radius = radius,
                             X = X, Y = Y, bin = bin,
                             bin.X = bin.X, bin.Y = bin.Y,
                             min.cells.per.region = min.cells.per.region)
    if(is.null(weights)) return(NULL)
    metacell <- tmpgcm %*% weights
    rownames(metacell) <- rownames(normdata)
    colnames(metacell) <- paste0(colnames(weights), "..", ct)
    return(metacell)
  }, mc.cores = ncores)
  metacell <- do.call(cbind, metacell_list)
  return(metacell)
}


GetKnnWeights <- function(metadata, k = 20, radius = 50,
                          X = "X", Y = "Y",
                          bin = "SpotID",
                          bin.X = "Spot.X",
                          bin.Y = "Spot.Y",
                          min.cells.per.region = 1){
  set.seed(39)
  require("Matrix")
  require("RANN")
  require("dplyr")
  locs <- as.data.frame(metadata[, c(X, Y, bin, bin.X, bin.Y)])
  colnames(locs) <- c("X", "Y", "SpotID", "Spot.X", "Spot.Y")
  spot_locs <- locs %>% distinct(SpotID, .keep_all = TRUE) %>% as.data.frame
  rownames(spot_locs) <- spot_locs$SpotID
  if(nrow(spot_locs)<3) return(NULL)

  ## sparse distance matrix
  if(k>nrow(locs)) k <- nrow(locs)-1
  knn <- RANN::nn2(data = as.matrix(locs[, 1:2]), query = as.matrix(spot_locs[,4:5]), k = k)
  weights <- as(matrix(0, nrow = nrow(locs), ncol = nrow(spot_locs)), "sparseMatrix")
  rownames(weights) <- rownames(locs)
  colnames(weights) <- rownames(spot_locs)

  if(length(weights)>1e8){
    ## Too big integers will be converted to NA
    for(i in 1:ncol(weights)){
      idx <- knn$nn.idx[i, ]
      idx <- idx[idx > 0]
      tmp <- knn$nn.dists[i, ]+1
      # tmp <- tmp[!(is.na(tmp)|is.infinite(tmp))]
      weights[idx, i] <- tmp
    }
    idx <- Matrix::which(weights>(radius+1), arr.ind = TRUE)
    for(j in unique(idx[,2])){
      i <- idx[idx[,2]==j, 1]
      weights[i, j] <- 0
    }
  }else{
    for(i in 1:k){
      idx <- nrow(weights) * (1:ncol(weights)-1) + knn$nn.idx[, i]
      idx <- idx[knn$nn.idx[, i]>0]
      tmp <- knn$nn.dists[, i]+1
      tmp <- tmp[!(is.na(tmp)|is.infinite(tmp))]
      weights[idx] <- tmp
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
