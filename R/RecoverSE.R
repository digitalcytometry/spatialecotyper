#' Recovery of SEs Using Pretrained NMF Models
#'
#' This function can recover SEs from scRNA-seq data or single-cell spatial data by assigning each single cell to an SE or NonSE.
#'
#' @param dat A numeric (sparse) matrix of gene expression data, from single-cell spatial
#' transcriptomics or scRNA-seq data.
#' @param celltypes Character vector specifying the cell type annotations for cells included in the gene expression `dat`.
#' If you're using the default model, cell types like B, CD4T, CD8T, NK, Plasma, Macrophage, DC, Fibroblast, Endothelial are expected.
#' @param scale Logical indicating whether to perform unit-variance normalization (default: TRUE). Change it with caution.
#' @param ncell.per.run Integer specifying the maximum number of cells per NMF prediction run to avoid memory issues.
#' @param Ws A list of cell-type-specific W matrices used to recover SE-specific
#' cell states. Each element in the list should be named after the corresponding
#' cell type and contain a W matrix from an NMF model.
#' @param ncores Integer specifying the number of CPU cores to use for parallel processing (default: 1).
#' @param se_results A list including a seurat object and a metadata with spatial
#' cluster annotations (SE column) returned by SpatialEcoTyper. When supplied, the `dat`
#' should be single cell gene expression data used for the SpatialEcoTyper analysis.
## @param .return.raw.prediction Default is FALSE. If TRUE, returns raw prediction scores instead of spatial ecotype assignments.
## Intended for advanced users.
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
                      ncores = 1,
                      se_results = NULL){
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
    assignRes$SE[assignRes$PredScore<0.6] <- "NonSE"
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
  # if(!is.null(se_results)){ ## SE recovery for single cell spatial data
  #   scmeta = se_results$metadata
  #   scmeta = scmeta[match(colnames(dat), rownames(scmeta)), ]
  #   rownames(scmeta) = colnames(dat)
  #   scmeta$State = paste0(scmeta$SE, "_", scmeta$CellType)
  #
  #   cts <- unique(intersect(names(Ws), unique(scmeta$CellType)))
  #   resDF <- lapply(cts, function(x){
  #     idx = (scmeta$CellType==x) & (!is.na(scmeta$SE))
  #     if(sum(idx)<2) return(NULL)
  #     tmpdat = dat[, idx]
  #     tmpmeta = scmeta[idx, ]
  #     if(scale) tmpdat = Seurat::ScaleData(tmpdat, verbose = FALSE)
  #     tmpdat = as.matrix(tmpdat)
  #     cell2cluster = matrix(0, nrow = ncol(tmpdat), ncol = length(unique(tmpmeta$SE)),
  #                           dimnames = list(colnames(tmpdat), unique(tmpmeta$SE)))
  #     idx = cbind(match(colnames(tmpdat), rownames(cell2cluster)),
  #                 match(tmpmeta$SE, colnames(cell2cluster)))
  #     idx = as.array(idx)
  #     cell2cluster[idx] = 1
  #     cell2cluster = t(t(cell2cluster) / colSums(cell2cluster))
  #     sedat = tmpdat %*% cell2cluster
  #
  #     H = NMFpredict(Ws[[x]], sedat, ncell.per.run = ncell.per.run,
  #                    scale = FALSE, ncores = ncores)
  #     preds = rownames(H)[apply(H, 2, which.max)]
  #     names(preds) = colnames(H)
  #     resDF = data.frame(Cluster = names(preds), SE = preds,
  #                         CellType = x, PredScore = apply(H, 2, max))
  #     return(resDF)
  #   })
  #   resDF <- do.call(rbind, resDF)
  #   resDF$State <- paste0(resDF$SE, "_", resDF$CellType)
  #   resDF$SE[resDF$PredScore<0.6] = "nonSE"
  #   if(flag=="default"){
  #     states <- readRDS(file.path(system.file("extdata", package = "SpatialEcoTyper"),
  #                                 "SE_CellStates.rds"))
  #     states <- setdiff(states, "SE01_B")
  #     resDF$State <- paste0(resDF$SE, "_", resDF$CellType)
  #     resDF$SE[!resDF$State %in% states] <- "nonSE"
  #     semap = c(SE01 = "SE1", SE02 = "SE2", SE03 = "SE3", SE04 = "SE4", SE05 = "SE5",
  #               SE06 = "nonSE", SE07 = "SE6", SE08 = "SE7", SE09 = "nonSE",
  #               SE10 = "SE8", SE11 = "SE9", nonSE = "nonSE")
  #     resDF$SE = semap[resDF$SE]
  #     resDF$State <- paste0(resDF$SE, "_", resDF$CellType)
  #   }
  #   preds <- resDF %>% count(Cluster, SE) %>%
  #     group_by(Cluster) %>% mutate(Frac = n / sum(n))
  #   wts <- 1 / table(preds$SE[preds$Cluster %in% scmeta$SE])
  #   preds$Frac <- preds$Frac * wts[preds$SE]
  #   preds <- preds %>% arrange(Cluster, -Frac) %>% distinct(Cluster, .keep_all = TRUE)
  #   idx <- match(scmeta$SE, preds$Cluster)
  #   preds <- preds$SE[idx]
  #   preds[is.na(preds)] = "nonSE"
  #   names(preds) <- rownames(scmeta)
  #   return(preds)
  # }
  # else{ ## SE recovery for Visium data
  #   statepreds <- mclapply(names(Ws), function(x){
  #     tmp = NMFpredict(W = Ws[[x]], dat, scale = scale)
  #     rownames(tmp) <- paste0(rownames(tmp), "_", x)
  #     tmp
  #   })
  #   statepreds <- do.call(rbind, statepreds)
  #   ses <- gsub("_.*", "", rownames(statepreds))
  #   state2se <- matrix(0, nrow = nrow(statepreds), ncol = length(unique(ses)),
  #                      dimnames = list(rownames(statepreds), unique(ses)))
  #   idx <- cbind(match(rownames(statepreds), rownames(state2se)),
  #                match(ses, colnames(state2se)))
  #   state2se[as.array(idx)] <- 1
  #   state2se <- t(t(state2se) / colSums(state2se))
  #   preds <- t(statepreds) %*% state2se
  #   preds <- preds / rowSums(preds)
  #
  #   if(flag=="default"){
  #     semap = c(SE01 = "SE1", SE02 = "SE2", SE03 = "SE3", SE04 = "SE4", SE05 = "SE5",
  #               SE06 = "NonSE", SE07 = "SE6", SE08 = "SE7", SE09 = "NonSE",
  #               SE10 = "SE8", SE11 = "SE9", NonSE = "NonSE")
  #     NonSE = rowSums(preds[, colnames(preds)%in%c("SE06", "SE09", "NonSE")])
  #     tmp = preds[, !colnames(preds)%in%c("SE06", "SE09", "NonSE")]
  #     colnames(tmp) = semap[colnames(tmp)]
  #     preds = cbind(tmp, NonSE = NonSE)
  #   }
  #   return(preds)
  # }
}
