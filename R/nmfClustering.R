#' Robust Clustering via NMF (non-negative matrix factorization)
#'
#' When one rank is provided, NMF clustering will be performed. When multiple ranks
#' are provided, this function will determine the optimal number of communities (K)
#' by assessing the cophenetic coefficient across different values of K.
#'
#' @param mat Numeric matrix (feature by sample) for clustering analysis.
#' @param ranks Numeric vector specifying the number of clusters to evaluate.
#' @param nrun.per.rank Integer specifying the number of runs per rank for clustering.
#' @param min.coph Numeric specifying the minimum cophenetic coefficient required for a rank to be optimal.
#' @param nmf.method Character string specifying the method for NMF analysis.
#' @param ncores Integer specifying the number of CPU cores to use for parallel processing.
#' @param plot Logical indicating whether to plot the results.
#' @param seed An integer used to seed the random number generator for NMF analysis.
#' @param ... Additional arguments to be passed to the nmf function.
#'
#' @return An NMFfitX1 object when only one rank is provided. A list containing
#' the optimal number of communities (bestK), a list of NMFfitX1 objects (NMFfits), and a ggplot
#' object (p) displaying the cophenetic coefficient across different values of K.
#'
#' @examples
#' library(NMF)
#' library(SpatialEcoTyper)
#' mat <- matrix(rnorm(1000, 3), 20)
#' mat[mat<0] = 0
#'
#' ## Specify one rank
#' result <- nmfClustering(mat = mat, ranks = 3, nrun.per.rank = 3)
#' predict(result)
#'
#' ## Determine optimal ranks by testing multiple ranks
#' result <- nmfClustering(mat = mat, ranks = 2:5, nrun.per.rank = 3)
#' result$p
#' result$bestK
#' predict(result$NMFfits[[paste0("K.", result$bestK)]])
#'
#' @export
#'
#' @import NMF
#' @importFrom parallel mclapply
#' @importFrom ggplot2 ggplot geom_point geom_line geom_vline theme theme_bw
#'
nmfClustering <- function(mat, ranks = 10,
                          nrun.per.rank = 30,
                          min.coph = 0.95,
                          nmf.method = "brunet",
                          ncores = 1, plot = FALSE,
                          seed = 2024, ...){
  ## Preprocess the data
  mat = as.matrix(mat)
  if(min(mat)<0) mat = posneg(mat)
  mat[is.na(mat)] = 0
  idx <- duplicated(rownames(mat))
  if(sum(idx)>0){
    rownames(mat)[!idx] <- paste0(rownames(mat)[!idx], "__pos")
    rownames(mat)[idx] <- paste0(rownames(mat)[idx], "__neg")
  }
  mat = mat[apply(mat, 1, var) > 0, ]

  ## Generate seeds for NMF analysis
  set.seed(seed)
  seeds = sample(1:2024, max(ranks)*nrun.per.rank, replace = TRUE)

  ## NMF analysis
  ranks = sort(ranks)
  if(length(ranks)>=nrun.per.rank){
    estim.list <- mclapply(ranks, function(k){
      tmp = lapply(1:nrun.per.rank, function(i){
        NMF::nmf(mat, k, nrun = 1, method = nmf.method, seed = seeds[k*i], ...)
      })
      suppressWarnings(NMF:::NMFfitX(tmp, .merge = T))
    }, mc.cores = ncores)
    names(estim.list) <- paste0("K.", ranks)
  }else{
    estim.list <- lapply(ranks, function(k){
      tmp = mclapply(1:nrun.per.rank, function(i){
        NMF::nmf(mat, k, nrun = 1, method = nmf.method, seed = seeds[k*i], ...)
      }, mc.cores = ncores)
      suppressWarnings(NMF:::NMFfitX(tmp, .merge = T))
    })
    names(estim.list) <- paste0("K.", ranks)
  }

  if(length(ranks)==1) return(estim.list[[1]])

  ## Post analysis
  cophs = data.frame(K = ranks, Cophenetic = unlist(lapply(estim.list, cophcor)))
  bestK = cophs$K[which(cophs$Cophenetic>min.coph)]
  if(length(bestK)<1){
    bestK = ranks[which.max(cophs$Cophenetic)]
  }else{
    bestK = max(bestK)
  }
  p <- ggplot(cophs, aes(x = K, y = Cophenetic)) +
    geom_point() + geom_line() +
    # geom_vline(xintercept = bestK, color = "red", linetype = "dashed") +
    theme(aspect.ratio = 1) +
    labs(x = "Number of communities", y = "Cophenetic coefficient") +
    theme_bw(base_size = 14)
  if(plot) plot(p)
  return(list(bestK = bestK, NMFfits = estim.list, p = p))
}
