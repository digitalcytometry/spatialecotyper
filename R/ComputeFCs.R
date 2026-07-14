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
#'
ComputeFCs <- function(normdata, scmeta, cluster = "SE",
                       Region = NULL, scale = FALSE,
                       ncores = 4){
  if(!"CellType" %in% colnames(scmeta)){
    stop("Metadata must include a column named 'CellType' for cell type annotations.")
  }
  if(!all(rownames(scmeta) %in% colnames(normdata))){
    stop("Rownames of `scmeta` do not match column names of `normdata`.")
  }
  scmeta$SE = scmeta[, cluster]
  scmeta <- scmeta[!is.na(scmeta$SE), ]
  normdata <- normdata[, match(rownames(scmeta), colnames(normdata)), drop = FALSE]

  ##### Within sample normalization ######
  lfcs <- mclapply(unique(scmeta$CellType), function(ct){
    if(sum(scmeta$CellType==ct)<10) return(NULL)
    tmpmeta = scmeta[scmeta$CellType==ct, ]
    tmpdata = normdata[, scmeta$CellType==ct, drop = FALSE]
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
    lfcs <- lfcs[, match(unique(scmeta$SE), colnames(lfcs)), drop = FALSE]
    colnames(lfcs) <- unique(scmeta$SE)
    lfcs[is.na(lfcs)] <- 0
    lfcs
  }, mc.cores = ncores)
  lfcs <- do.call(rbind, lfcs)
  return(lfcs)
}
