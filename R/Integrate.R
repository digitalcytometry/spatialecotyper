#' Integrate Spatial Clusters From Multiple Samples Via Similarity Network Fusion
#'
#' This function identifies conserved spatial ecotypes based on cell type specific
#' gene expression signature of spatial clusters from all samples.
#'
#' @param avgexprs Gene expression signature of spatial clusters, where each column represents a spatial cluster.
#' @param Region A character vector specifying the region annotations for all spatial clusters.
#' @param nfeatures Integer specifying the maximum number of top variable genes to select for each cell type.
#' @param min.features Integer specifying the minimum number of shared features (genes) required across samples.
#' @param minibatch Integer specifying the number of columns to process in each minibatch in the SNF analysis.
#' Default is 5000. This option splits the matrix into smaller chunks (minibatch), thus reducing memory usage.
#' @param ncores Integer specifying the number of cores for parallel processing. Default is 1.
#'
#' @return Integrated similarity matrix of spatial clusters across all samples.
#'
Integrate <- function(avgexprs, Region = NULL,
                      nfeatures = 3000,
                      min.features = 5,
                      minibatch = 5000,
                      ncores = 1){
  #### Cell-type specific similarity matrix ####
  celltypes <- gsub("\\.\\..*", "", rownames(avgexprs))
  if(!is.null(Region)) names(Region) <- colnames(avgexprs)
  cors_rank = mclapply(unique(celltypes), function(ct){
    tmpdat = avgexprs[celltypes==ct, ]
    tmpdat = tmpdat[, colSums(!is.na(tmpdat))>min.features]
    tmpdat = tmpdat[, colSums(tmpdat!=0, na.rm = TRUE)>min.features]
    # Select variable genes
    samples = gsub("\\.\\..*", "", colnames(tmpdat))
    if(length(table(samples)) < 2 | min(table(samples))<3) return(NULL)
    vars = mclapply(unique(samples), function(ss){
      apply(tmpdat[, samples==ss], 1, var, na.rm = TRUE)
    })
    vars = do.call(cbind, vars)
    vars = vars[apply(vars, 1, min)>0, ]
    var.ranks = apply(-vars, 2, rank)
    vars = rownames(var.ranks)[order(rowMeans(log(var.ranks)))]
    vars = vars[1:min(length(vars), nfeatures)]
    if(length(vars) < min.features) return(NULL)
    tmpdat = tmpdat[vars, ]
    tmpcor = suppressWarnings(cor(as.matrix(tmpdat), method = "spearman"))
    tmpcor[is.na(tmpcor)] = 0
    tmpcor[tmpcor<0] = 0

    if(!is.null(Region)){ # Add weights to spatial clusters in the same region
      tmpregion = Region[rownames(tmpcor)]
      for(g in unique(tmpregion)){
        idx = tmpregion==g
        tmpcor[idx, idx] = tmpcor[idx, idx]*1.1
      }
    }

    samples = gsub("\\.\\..*", "", rownames(tmpcor))
    tmpcor = mclapply(unique(samples), function(ss){
      sscor = apply(tmpcor[samples==ss, ], 2, rank)
      sscor = t(t(sscor) / colSums(sscor))
      sscor
    }, mc.cores = ncores)
    tmpcor = do.call(rbind, tmpcor)
    tmpcor = tmpcor[, rownames(tmpcor)]
    tmpcor = tmpcor + t(tmpcor)
    return(tmpcor)
  }, mc.cores = ncores)
  names(cors_rank) <- unique(celltypes)
  cors_rank = cors_rank[lengths(cors_rank)>0]
  cors_rank = fillspots(cors_rank)
  ## SNF of correlation matrices
  obj = SNF2(cors_rank, ncores = ncores,
             minibatch = minibatch, verbose = TRUE)
  obj = rankSparse(obj)
  obj = (obj + t(obj)) / 2
  return(obj)
}
