#' Recovery of SEs Using Pretrained NMF Models
#'
#' This function can recover SEs from scRNA-seq data or single-cell spatial data
#' by assigning each single cell to an SE or NonSE.
#'
#' @param dat A numeric (sparse) matrix of gene expression data, from single-cell spatial
#' transcriptomics or scRNA-seq data.
#' @param celltypes Character vector specifying the cell type annotations for
#' cells included in the gene expression `dat`. If you're using the default model,
#' cell types including B, CD4T, CD8T, NK, Plasma, Macrophage, DC, Fibroblast, Endothelial are expected.
#' @param scale Logical indicating whether to perform unit-variance normalization
#' (default: TRUE). Change it with caution.
#' @param ncell.per.run Integer specifying the maximum number of cells per NMF
#' prediction run to avoid memory issues.
#' @param Ws A list of cell-type-specific W matrices used to recover SE-specific
#' cell states. Each element in the list should be named after the corresponding
#' cell type and contain a W matrix from an NMF model.
#' @param min.score A numeric threshold (0-1) specifying the minimum prediction
#' score for SE classification; cells with lower scores are assigned to NonSE.
#' @param ncores Integer specifying the number of CPU cores to use for parallel processing (default: 1).
#'
#'
#' @return A data frame with four columns, including CID (cell id), CellType (cell type),
#' SE (spatial ecotype), and PredScore (prediction probability for the assigned SE).
#'
#' @examples
#' # see https://digitalcytometry.github.io/spatialecotyper/articles/Recovery_scRNA.html
#' # see https://digitalcytometry.github.io/spatialecotyper/articles/Recovery_Spatial.html
#'
#' @export
#'
#' @importFrom parallel mclapply
#'

RecoverSE <- function(dat, celltypes = NULL,
                      scale = TRUE,
                      ncell.per.run = 500,
                      Ws = NULL,
                      min.score = 0.6,
                      ncores = 1){
  ## deprecated: .return.raw.prediction = FALSE
  flag = ifelse(is.null(Ws), "default", "custom")
  ## Load NMF models
  if(is.null(Ws)){
    Wfiles <- list.files(system.file("extdata", package = "SpatialEcoTyper"),
                         "MERSCOPE_W_v1.*.rds", full.names = TRUE)
    Ws <- lapply(Wfiles, readRDS)
    names(Ws) <- gsub(".*v1_|.rds", "", Wfiles)
  }
  if(!scale) warning("Unit-variance normalization is essential for the prediction.")
  if(!is.null(celltypes)){ ## SE recovery for scRNA-seq data
    if(ncol(dat)!=length(celltypes)){
      stop("The cell type annotations do not match the columns of expression data")
    }
    if(!is.null(names(celltypes))){
      celltypes = celltypes[names(celltypes) %in% colnames(dat)]
      dat = dat[, colnames(dat) %in% names(celltypes)]
      celltypes = celltypes[match(colnames(dat), names(celltypes))]
    }
    names(celltypes) = colnames(dat)
    cts <- unique(intersect(names(Ws), unique(celltypes)))
    ses = sort(unique(unlist(lapply(Ws, colnames))))
    if(length(cts)==0){
      stop("At least one of ", paste0(names(Ws), collapse = ", "),
           " should be included for SE recovery")
    }
    if(all(dat>=0) & max(dat)>80){
      message("automated log2 normalization")
      dat = log2(dat+1)
    }
    resDF <- lapply(cts, function(x){
      message(Sys.time(), " Recover SE cell states from ", x)
      if(sum(celltypes==x)<20) return(NULL)
      tmpdat <- dat[, celltypes==x]
      H <- NMFpredict(Ws[[x]], tmpdat, ncell.per.run = ncell.per.run,
                      scale = scale, ncores = ncores, sum2one = TRUE)
      H = H[, match(ses, colnames(H))]
      colnames(H) = ses
      return(H)
    })
    resDF <- do.call(rbind, resDF)

    idx <- !(colnames(dat) %in% rownames(resDF))
    if(sum(idx)>0){## add unpredicted cells
      tmp = matrix(0, nrow = sum(idx), ncol = ncol(resDF))
      rownames(tmp) = colnames(dat)[idx]
      colnames(tmp) = colnames(resDF)
      resDF = rbind(resDF, tmp)
    }
    resDF[is.na(resDF)] = 0

    # if(.return.raw.prediction){
    #   return(resDF)
    # }
    assignRes = data.frame(CID = rownames(resDF),
                           CellType = celltypes[match(rownames(resDF), names(celltypes))],
                           SE = colnames(resDF)[apply(resDF, 1, which.max)],
                           PredScore = apply(resDF, 1, max))
    assignRes$SE[assignRes$PredScore<min.score] <- "NonSE"
    if(flag=="default"){
      states <- readRDS(file.path(system.file("extdata", package = "SpatialEcoTyper"),
                                  "SE_CellStates.rds"))
      tmpstat <- paste0(assignRes$SE, "_", assignRes$CellType)
      assignRes$SE[!tmpstat %in% states] <- "NonSE"
      semap = c(SE01 = "SE1", SE02 = "SE2", SE03 = "SE3", SE04 = "SE4", SE05 = "SE5",
                SE06 = "NonSE", SE07 = "SE6", SE08 = "SE7", SE09 = "NonSE",
                SE10 = "SE8", SE11 = "SE9", NonSE = "NonSE")
      assignRes$SE = semap[assignRes$SE]
      assignRes <- assignRes[match(colnames(dat), assignRes$CID), ]
      assignRes$SE[is.na(assignRes$SE)] = "NonSE"
    }
    return(assignRes)
  }
}
