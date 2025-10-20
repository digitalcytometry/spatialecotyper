#' Generate Principal Component (PC) List for Spatial Neighborhoods
#'
#' This function performs principal component analysis based on cell-type-specific
#' gene expression of spatial neighborhoods.
#'
#' @param mergedncem A matrix of cell-type-specific gene expression data, with rows
#' representing genes and columns representing spatial neighborhoods by cell type.
#' @param min.cells Integer, minimum number of non-zero counts required per gene to retain the gene in the analysis. Default is 3.
#' @param min.features Integer, minimum number of features (genes) required per
#' spatial neighborhood to retain the neighborhood in the analysis. Default is 5.
#' @param nfeatures Integer, number of variable features to select for PCA. Default is 3000.
#' @param ncores Integer, number of cores to use for parallel processing. Default is 1.
#' @param do.scale Logical, whether to scale data before performing PCA. Default is TRUE.
#'
#' @details
#' The function first filters the input matrix based on the specified minimum number
#' of cells and features. It then creates Seurat objects for each cell type, normalizes
#' the data, identifies variable features, scales the data, and performs PCA.
#' The resulting PCs for each cell type are returned as a list.
#'
#' @return A named list of matrices, where each matrix contains the PCs for a specific
#' cell type. Cell types with insufficient data are excluded from the result.
#'
#' @import parallel Seurat
#' @export
#'
GetPCList <- function(mergedncem,
                      min.cells = 3,
                      min.features = 5,
                      nfeatures = 3000,
                      ncores = 1, do.scale = TRUE){
  celltypes <- gsub(".*\\.+", "", colnames(mergedncem))
  emb_list <- parallel::mclapply(unique(celltypes), function(ct){
    ncem = mergedncem[, which(celltypes==ct), drop = FALSE]
    colnames(ncem) = gsub("\\.+.*", "", colnames(ncem))
    if(sum(rowSums(ncem>0)>=min.cells)<min.features) return(NULL)
    ncem = ncem[rowSums(ncem>0)>=min.cells, ]
    if(sum(colSums(ncem>0)>=min.features)<max(min.cells, 5)) return(NULL)
    ncem = ncem[, colSums(ncem>0)>=min.features]
    tmpobj = CreateSeuratObject(ncem - rowMeans(ncem))
    seurat_version = as.integer(gsub("\\..*", "", as.character(packageVersion("Seurat"))))
    if(seurat_version>=5) tmpobj[["RNA"]]$data = tmpobj[["RNA"]]$counts
    tmpobj = tmpobj %>%
      FindVariableFeatures(selection.method = "dispersion",
                           nfeatures = nfeatures, verbose = FALSE) %>%
      ScaleData(do.scale = do.scale, verbose = FALSE) %>%
      RunPCA(verbose = FALSE, npcs = min(ncol(ncem)-3, 30))
    t(Embeddings(tmpobj, "pca"))
  }, mc.cores = ncores)
  names(emb_list) = unique(celltypes)

  idx <- unlist(lapply(emb_list, is.null))
  if(sum(idx)>0){
    warning("\t\tExclude ", paste0(unique(celltypes)[idx], collapse = ", "),
            " due to limited number of spatial metacells\n")
    emb_list <- emb_list[!idx]
  }
  emb_list
}
