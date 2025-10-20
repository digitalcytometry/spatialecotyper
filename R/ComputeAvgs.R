#' Compute Cell-Type-Specific Average Expression of Spatial Clusters
#'
#' @param normdata Numeric matrix of normalized expression data, where rows
#' represent genes and columns represent cells.
#' @param scmeta Data frame containing metadata associated with each cell,
#' including spatial cluster and cell type annotations.
#' @param cluster Character string specifying the column name in 'scmeta'
#' containing spatial cluster annotations.
#' @param Region Character string specifying the column name in metadata data
#' frames containing region annotations (default: NULL).
#' @param scale A boolean specifying whether to do univariance normalization.
#' @param ncores Integer specifying the number of CPU cores to use for parallel processing.
#'
#' @return A matrix of average expression, where rows represent genes and columns
#' represent spatial clusters from the sample.
#'
#' @examples
#' library(data.table)
#' library(Seurat)
#' library(SpatialEcoTyper)
#' library(googledrive)
#' drive_deauth() # no Google sign-in is required
#' drive_download(as_id("1CoQmU3u8MoVC8RbLUvTDQmOuJJ703HHB"),
#'               "HumanMelanomaPatient1_subset_counts.tsv.gz", overwrite = TRUE)
#' drive_download(as_id("1nSPj2zRywFUdbo1fwiz77ds4NuM6bmV2"),
#'               "Melanoma1_subset_SpatialEcoTyper_results.rds", overwrite = TRUE)
#' scdata <- fread("HumanMelanomaPatient1_subset_counts.tsv.gz",
#'                 sep = "\t",header = TRUE, data.table = FALSE)
#' rownames(scdata) <- scdata[, 1]
#' scdata <- as.matrix(scdata[, -1])
#' tmpobj <- CreateSeuratObject(scdata) %>%
#'         SCTransform(clip.range = c(-10, 10), verbose = FALSE)
#' seurat_version = as.integer(gsub("\\..*", "", as.character(packageVersion("SeuratObject"))))
#' if(seurat_version<5){
#'   normdata <- GetAssayData(tmpobj, "data")
#' }else{
#'   normdata <- tmpobj[["SCT"]]$data
#' }
#' metadata = readRDS("Melanoma1_subset_SpatialEcoTyper_results.rds")$metadata
#'
#' # Construct cell-type-specific gene expression signatures of SEs
#' avgexprs <- ComputeAvgs(normdata = normdata, scmeta = metadata)
#' head(avgexprs)
#'
#' @import Seurat
#' @keywords internal
#'
ComputeAvgs <- function(normdata, scmeta, cluster = "SE",
                       Region = NULL, scale = TRUE, ncores = 4){
  scmeta$SE = scmeta[, cluster]
  scmeta = scmeta[!is.na(scmeta$SE), ]
  cids = na.omit(intersect(rownames(scmeta), colnames(normdata)))
  scmeta = scmeta[cids, ]
  normdata = normdata[, cids]
  if(max(normdata)>50 & min(normdata)>=0) normdata = log1p(normdata+1)
  ##### Within sample normalization ######
  lfcs <- mclapply(unique(scmeta$CellType), function(ct){
    if(sum(scmeta$CellType==ct)<10) return(NULL)
    tmpdata = normdata[, scmeta$CellType==ct]
    tmpmeta = scmeta[scmeta$CellType==ct, ]
    if(length(unique(tmpmeta$SE))<2) return(NULL)
    if(scale){
      if(!is.null(Region) && (Region%in%colnames(scmeta))){
        if(seurat_version>=5){
          tmpdat <- Znorm(tmpdata, groups = tmpmeta[, Region])
        }else{
          tmpdat <- Znorm(tmpdata, groups = tmpmeta[, Region])
        }
      }else{
        tmpdat <- ScaleData(tmpdata, verbose = FALSE)
      }
    }
    tmp = Matrix(0, nrow = ncol(tmpdata),
                 ncol = length(unique(tmpmeta$SE)), sparse = TRUE)
    colnames(tmp) = sort(unique(tmpmeta$SE))
    idx = cbind(1:nrow(tmpmeta), match(tmpmeta$SE, colnames(tmp)))
    tmp[as.array(idx)] = 1
    tmp = Matrix::t(Matrix::t(tmp) / Matrix::colSums(tmp))
    lfcs <- tmpdata %*% tmp
    rownames(lfcs) <- paste0(ct, "..", rownames(lfcs))
    lfcs
  }, mc.cores = ncores)
  names(lfcs) <- unique(scmeta$CellType)
  lfcs = lfcs[lengths(lfcs)>0]
  lfcs <- mclapply(lfcs, function(x){
    x = as.matrix(x)
    x = x[, match(unique(scmeta$SE), colnames(x))]
    colnames(x) = unique(scmeta$SE)
    x
  }, mc.cores = ncores)
  lfcs <- do.call(rbind, lfcs)
  return(lfcs)
}


#' Compute Cell-Type-Specific Fold Changes (FCs) for Spatial Clusters
#'
#' This function computes fold changes (FCs) for spatial transcriptomic data,
#' comparing expression levels between spatial clusters within each cell type.
#'
#' @param normdata Numeric matrix of normalized expression data, where rows
#' represent genes and columns represent cells.
#' @param scmeta Data frame containing metadata associated with each cell,
#' including spatial cluster and cell type annotations.
#' @param cluster Character string specifying the column name in 'scmeta'
#' containing spatial cluster annotations.
#' @param Region Character string specifying the column name in metadata data
#' frames containing region annotations (default: NULL).
#' @param scale A boolean specifying whether to do univariance normalization.
#' @param ncores Integer specifying the number of CPU cores to use for parallel processing.
#'
#' @return A matrix of fold changes (FCs), where rows represent genes and columns
#' represent spatial clusters from the sample.
#'
#'
#' @import Seurat
#' @export
#'
ComputeFCs <- function(normdata, scmeta, cluster = "SE",
                       Region = NULL, scale = FALSE,
                       ncores = 4){
  if(!"CellType" %in% colnames(scmeta)){
    stop("the meta data have to include a column (CellType) for cell type annotations")
  }
  scmeta$SE = scmeta[, cluster]
  scmeta <- scmeta[!is.na(scmeta$SE), ]
  normdata <- normdata[, match(rownames(scmeta), colnames(normdata))]

  ##### Within sample normalization ######
  lfcs <- mclapply(unique(scmeta$CellType), function(ct){
    if(sum(scmeta$CellType==ct)<10) return(NULL)
    tmpmeta = scmeta[scmeta$CellType==ct, ]
    tmpdata = normdata[, scmeta$CellType==ct]
    if(scale){
      if(!is.null(Region) && (Region%in%colnames(scmeta))){
        tmpdata <- Znorm(tmpdata, groups = tmpmeta[, Region])
      }else{
        tmpdata <- ScaleData(tmpdata)
      }
    }
    if(length(unique(tmpmeta$SE))<2) return(NULL)
    lfcs <- lapply(unique(tmpmeta$SE), function(cl){
      rowMeans(tmpdata[, tmpmeta$SE==cl, drop=FALSE]) - rowMeans(tmpdata[, tmpmeta$SE!=cl, drop=FALSE])
    })
    lfcs <- do.call(cbind, lfcs)
    colnames(lfcs) <- unique(tmpmeta$SE)
    rownames(lfcs) <- paste0(ct, "..", rownames(lfcs))
    lfcs <- lfcs[, match(unique(scmeta$SE), colnames(lfcs))]
    colnames(lfcs) <- unique(scmeta$SE)
    lfcs[is.na(lfcs)] <- 0
    lfcs
  }, mc.cores = ncores)
  lfcs <- do.call(rbind, lfcs)
  return(lfcs)
}
