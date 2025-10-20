#' Enhanced Similarity Network Fusion
#'
#' Similarity Network Fusion (SNF) integrates multiple views (similarity matrices) to
#' construct an overall status matrix. This function is adopted from the SNFtool
#' (https://github.com/maxconway/SNFtool/) and has been enhanced for unsupervised analysis of spatial ecosystem.
#' The new function supports sparse matrix and missing data in the input matrices.
#'
#' @param Wall List of similarity matrices. Each element of the list is a square,
#' symmetric matrix that shows affinities of the data points from a certain view.
#' @param K Number of neighbors in K-nearest neighbors part of the algorithm.
#' @param t Number of iterations for the diffusion process.
#' @param minibatch Integer specifying the number of columns to process in each minibatch.
#' Default is 5000. This option splits the matrix into smaller chunks (minibatch),
#' thus reducing memory usage.
#' @param ncores Integer specifying the number of CPU cores to use for parallel processing.
#' @param verbose Boolean specifying whether to show the progress messages.
#' @return A fused matrix.
#'
#' @import Matrix
#' @export
#'
SNF2 <- function(Wall, K = 10, t = 10,
                 minibatch = 5000, ncores = 4,
                 verbose = FALSE){ #
  require("Matrix")
  # Similarity Network Fusion takes multiple views of a network (Wall) and
  # fuses them together to create a overall affinity matrix.
  #
  # Args:
  #   Wall: List of matrices, each element is a square symmetric affinity
  #       matrix.
  #   K: Number of neighbors used in the K-nearest neighbours step
  #   t: Number of iterations for the diffusion process
  #
  # Returns:
  #   W: Unified similarity graph of all data types in Wall.

  check_wall_names <- function(Wall){
    name_match <- function(names_A, names_B){
      return(identical(dimnames(names_A), dimnames(names_B)))
    }

    return(all(unlist(lapply(Wall, FUN=name_match, Wall[[1]]))))
  }

  #Check if Wall names are consistant across all matrices in Wall
  wall.name.check <- check_wall_names(Wall)
  wall.names <- dimnames(Wall[[1]])
  if(!wall.name.check){
    warning("Dim names not consistent across all matrices in Wall.
            Returned matrix will have no dim names.")
  }

  LW <- length(Wall)

  #Normalization method for affinity matrices
  normalize <- function(X){
    if(is(X, "sparseMatrix")){
      row.sum.mdiag <- Matrix::rowSums(X) - diag(X)
    }else{
      row.sum.mdiag <- rowSums(X) - diag(X)
    }
    #If rowSumx(X) == diag(X), set row.sum.mdiag to 1 to avoid div by zero
    row.sum.mdiag[row.sum.mdiag == 0] <- 1
    X <- X / (2*(row.sum.mdiag))
    diag(X) <- 0.5
    if(is(X, "sparseMatrix")){
      X <- (X + Matrix::t(X))/2
    }else{
      X <- (X + t(X))/2
    }
    return(X)
  }

  #Normalize different networks to avoid scale problems.
  if(verbose) message(Sys.time(), " Normalize networks ...")
  Wall <- mclapply(Wall, normalize, mc.cores = ncores)

  ### Calculate the local transition matrix.
  if(verbose) message(Sys.time(), " Calculate the local transition matrix ...")
  newW <- lapply(Wall, function(X){ (.dominateset(X, K, ncores)) })
  newW <- mclapply(newW, normalize, mc.cores = ncores)

  #Perform the diffusion for t iterations
  if(verbose) message(Sys.time(), " Perform the diffusion ...")

  # sum_matrix <- function(Mats){ ## Sum up a list of matrices
  #   sumres = Mats[[1]]
  #   for (i in 2:length(Mats)) {
  #     sumres <- sumres + Mats[[i]]
  #   }
  #   return(sumres)
  # }

  for (i in 1:t){
    if(verbose) message("\t", Sys.time(), " Iteration: ", i)
    Wall <- mclapply(1:LW, function(j){
      sumWJ <- Reduce("+", Wall[-j]) / (LW-1)
      x <- matrixMultiply(newW[[j]], sumWJ, minibatch = minibatch)
      x <- matrixMultiply(x, Matrix::t(newW[[j]]), minibatch = minibatch)
      normalize(x)
    }, mc.cores = ncores)
  }

  LW <- Reduce("+", mclapply(Wall, function(x) { x>0 }, mc.cores = ncores))
  LW <- as(LW, "sparseMatrix")

  Wall <- Reduce("+", Wall) / LW
  Wall@x[is.infinite(Wall@x)] <- 0
  Wall@x[is.nan(Wall@x)] <- 0
  rownames(Wall) = wall.names[[1]]
  colnames(Wall) = wall.names[[2]]
  return(Wall)
}

.dominateset <- function(xx, KK=20, ncores = 8){
  require("Matrix")
  if(nrow(xx)>5000) xx <- as(xx, "sparseMatrix")
  ###This function outputs the top KK neighbors.
  zero <- function(x) {
    s = sort(x, index.return=TRUE)
    x[s$ix[1:(length(x)-KK)]] = 0
    return(x)
  }
  # normalize <- function(X) X / Matrix::rowSums(X)
  xx <- parallel::mclapply(1:nrow(xx), function(i) {
    zero(xx[i,])
  }, mc.cores = ncores)
  xx <- do.call(rbind, xx)
  # xx <- normalize(xx)
  xx <- as(xx, "sparseMatrix")
  return(xx)
}

#' Matrix Multiplication with Minibatching and Parallel Processing
#'
#' This function performs matrix multiplication of two matrices (`mat1` and `mat2`)
#' using minibatching to manage memory usage and parallel processing to speed up
#' computation. It is particularly useful for large matrices where full multiplication
#' would otherwise be computationally intensive or memory prohibitive.
#'
#' @param mat1 A matrix with dimensions (m x n).
#' @param mat2 A matrix with dimensions (n x p).
#' @param minibatch The number of columns from `mat2` to process in each minibatch. Default is 5000.
#' @param ncores The number of cores to use for parallel processing. Default is 1 (no parallel processing).
#'
#' @return A matrix with dimensions (m x p) representing the product of `mat1` and `mat2`.
#' @examples
#' # Example usage:
#' mat1 <- matrix(runif(1000), nrow=100, ncol=10)
#' mat2 <- matrix(runif(2000), nrow=10, ncol=200)
#' result <- matrixMultiply(mat1, mat2, minibatch=100, ncores=2)
#' @import parallel
#' @export
matrixMultiply <- function(mat1, mat2, minibatch = 5000, ncores = 1){
  cycles = ceiling(ncol(mat2) / minibatch)
  res = mclapply(1:cycles, function(ii){
    query_idx = ((ii-1)*minibatch+1):min(ii*minibatch, ncol(mat2))
    x = mat1 %*% mat2[, query_idx, drop = FALSE]
    return(x)
  }, mc.cores = ncores)
  do.call(cbind, res)
}
