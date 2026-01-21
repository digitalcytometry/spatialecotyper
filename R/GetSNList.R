#' Construct Similarity Network
#'
#' This function generates a similarity network (SN) for a given embeddings.
#' Specifically, the function computes K-nearest neighbor (KNN) graphs based
#' on the provided embedding.
#'
#' @param emb A matrix contains the embedding, with rows correspond to principal
#' components (PCs) or features, and columns correspond to spatial neighborhoods.
#' @param npcs An integer specifying the number of principal components or features
#' for the analysis. Default is 20.
#' @param k An integer specifying the number of nearest neighbors to consider when constructing the KNN graph.
#' Default is 50.
#'
#' @return A matrix for the similarity network (SN). The rows and columns represent spatial neighborhood,
#' and the values represent inverted distance.
#'
#' @details
#' The function computes a KNN graph using the input embedding and inverts the resulting
#' distances to similarity scores.
#'
#' @importFrom RANN nn2
#' @importFrom Matrix Matrix drop0
#' @importFrom parallel mclapply
#' @export
getSN <- function(emb, npcs = 20, k = 50){
  tmpK = min(k, ncol(emb)-1)
  knn = RANN::nn2(data = t(emb), query = t(emb), k = tmpK)
  W = Matrix(0, nrow = ncol(emb), ncol = ncol(emb), sparse = TRUE)
  rownames(W) = colnames(emb)
  colnames(W) = colnames(emb)
  if(length(W)>1e8){
    ## Too big integers will be converted to NA
    for(i in 1:ncol(W)){
      idx = knn$nn.idx[i, ]
      idx = idx[idx > 0]
      tmp = knn$nn.dists[i, ]+1
      # tmp <- tmp[!(is.na(tmp)|is.infinite(tmp))]
      W[idx, i] = tmp
    }
  }else{
    for(i in 1:tmpK){
      idx = nrow(W) * (1:ncol(W)-1) + knn$nn.idx[, i]
      idx = idx[knn$nn.idx[, i]>0]
      tmp = knn$nn.dists[, i]+1
      tmp = tmp[!(is.na(tmp)|is.infinite(tmp))]
      W[idx] = tmp
    }
  }
  W = drop0(W)
  W@x = 1 / (W@x+1)
  return(W)
}


#' Construct Cell-Type-Specific Similarity Network
#'
#' This function generates a list of similarity networks (SNs) for a given list of embeddings.
#' The embeddings represent spatial neighborhoods, and the function computes K-nearest neighbor (KNN)
#' graphs based on the provided embeddings to create similarity matrices.
#'
#' @param emb_list A list of matrices where each matrix contains the embeddings for a single sample.
#' The rows correspond to different principal components (PCs) or features, and columns correspond to spatial regions or cells.
#' @param npcs An integer specifying the number of principal components or features to use from each embedding matrix.
#' Default is 20.
#' @param k An integer specifying the number of nearest neighbors to consider when constructing the KNN graph.
#' Default is 50.
#' @param min.cts.per.region An integer specifying the minimum number of distinct cell types that must be present
#' in each spatial region for it to be retained in the analysis. Default is 1.
#' @param ncores An integer specifying the number of cores to use for parallel computation. Default is 1.
#'
#' @return A list of sparse matrices where each matrix represents a similarity network (SN) for a sample.
#' The rows and columns of each matrix correspond to spatial regions, and the values represent similarity scores.
#'
#' @details
#' This function first trims the input embeddings to retain only the top `npcs` principal components or features.
#' It then filters out spatial regions that do not meet the minimum cell type threshold.
#' For each sample, the function computes a KNN graph using the trimmed embeddings and converts the resulting
#' distances to similarity scores. Missing spatial regions are filled with zeros to ensure all samples have the same regions.
#'
#' @importFrom RANN nn2
#' @importFrom Matrix Matrix drop0
#' @importFrom parallel mclapply
#' @export
GetSNList <- function(emb_list,
                      npcs = 20,
                      k = 50,
                      min.cts.per.region = 1,
                      ncores = 1){
  emb_list <- lapply(emb_list, function(x){
    x[1:min(npcs, nrow(x)), ]
  })
  nspots <- table(unlist(lapply(emb_list, colnames)))
  if(sum(nspots<min.cts.per.region)>0){
    message("\t\tRemove ", sum(nspots<min.cts.per.region),
            " spatial neighborhoods with fewer than ",
            min.cts.per.region, " cell types ")
    nspots <- names(nspots)[nspots>=min.cts.per.region]
    emb_list <- lapply(emb_list, function(x){
      x[, colnames(x) %in% nspots]
    })
  }
  snlist <- mclapply(emb_list, function(emb){
    getSN(emb, npcs = npcs, k = k)
  }, mc.cores = ncores)
  snlist = fillspots(snlist)
  return(snlist)
}

#' Handle Missing Values
#'
#' This function fills missing spatial regions in similarity network matrices with zeros.
#' It ensures that all matrices in the list have the same set of spatial regions.
#'
#' @param snlist A list of similarity network matrices where the rows and columns correspond to spatial regions.
#'
#' @return A list of similarity network matrices with missing spatial regions filled in with zero values.
#'
#' @details
#' The function first identifies all unique spatial regions across the input list of matrices.
#' It then ensures that each matrix has the same set of regions by adding missing regions
#' with zero similarity values. This step is crucial for downstream analyses that require consistent
#' dimensions across samples.
#'
#' @export
fillspots <- function(snlist){
  ### consider missing cell types as 0
  spots <- unique(unlist(lapply(snlist, colnames)))
  lapply(snlist, function(x){
    miss = setdiff(spots, colnames(x))
    if(length(miss)>0){
      tmp = matrix(0, nrow = nrow(x), ncol = length(miss))
      colnames(tmp) = miss
      x = cbind(x, tmp)
      tmp = matrix(0, nrow = length(miss), ncol = ncol(x))
      rownames(tmp) = miss
      x = rbind(x, tmp)
      x
    }
    x[match(spots, rownames(x)), match(spots, colnames(x))]
  })
}


# GetSNList <- function(emb_list,
#                       method = "euclidean",
#                       npcs = 20,
#                       k = 20,
#                       min.cts.per.region = 1,
#                       ncores = 1){
#   emb_list <- lapply(emb_list, function(x){
#     x[1:min(npcs, nrow(x)), ]
#   })
#   nspots <- table(unlist(lapply(emb_list, colnames)))
#   if(sum(nspots<min.cts.per.region)>0){
#     message("\t\tRemove ", sum(nspots<min.cts.per.region),
#             " spatial neighborhoods with fewer than ",
#             min.cts.per.region, " cell types ")
#     nspots <- names(nspots)[nspots>=min.cts.per.region]
#     emb_list <- lapply(emb_list, function(x){
#       x[, colnames(x) %in% nspots]
#     })
#   }
#   snlist <- mclapply(emb_list, function(x){
#     dists <- as.matrix(stats::dist(t(x), method = method))
#     ## Replace distance of 0 with a random small distance
#     # dists[dists==0] <- quantile(dists[dists>0], 0.01)
#     W <- affinityMatrix(dists, K = min(k, nrow(dists)-1), sigma = 0.5)
#     rownames(W) <- colnames(x)
#     colnames(W) <- colnames(x)
#     return(W)
#   }, mc.cores = ncores)
#   rm(emb_list)
#   snlist <- fillspots(snlist)
#   return(snlist)
# }
#
