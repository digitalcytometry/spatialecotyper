#' Train Cell Type-Specific NMF Models for Recovering Spatial EcoTypes
#'
#' This function trains cell type-specific NMF (Non-Negative Matrix Factorization)
#' models to recover SE-specific cell states from single-cell data, as part of the
#' Spatial EcoTyper analysis workflow.
#' It downsamples cells for training when the dataset size is large, and selects a
#' subset of features with the highest specificity.
#'
#' @param scdata Numeric matrix containing single-cell expression data.
#' @param scmeta Data frame containing metadata information associated with single-cell data,
#' including cell types and spatial clusters.
#' @param CellType Character string specifying the column name in the metadata data
#' frame containing cell type annotations.
#' @param SE Character string specifying the column name in the metadata data frame
#' containing spatial ecotype annotations.
#' @param scale Boolean specifying whether to perform univariance normalization
#' before training the models (default: TRUE).
#' @param Sample Character string specifying the column name in the metadata data frame
#' containing sample annotations. If specified, the univariance normalization will be
#' performed within each sample.
#' @param balance.sample Boolean specifying whether to perform balance the cells from all samples
#' before training the models (default: TRUE).
#' @param nfeature Integer specifying the top variable features for training the models (default: 2000).
#' @param nfeature.per.se Integer specifying the maximal number of features to select for each SE (default: 50).
#' @param min.cells Integer specifying the minimal number of cells required for
#' each SE cell state.
#' @param downsample Integer specifying the number of cells per cell type (downsampling)
#' for training the NMF models (default: 2500).
#' @param ncores Integer specifying the number of CPU cores to use for parallel processing.
#' @param seed Integer specifying the seed for random sampling during downsampling.
#'
#' @return A list of cell type-specific NMF models, each represented by its corresponding
#' factorization matrix W.
#'
#' @export
#'
#' @importFrom parallel mclapply
#'
NMFGenerateWList <- function(scdata, scmeta,
                             CellType = "CellType",
                             SE = "SE",
                             scale = TRUE,
                             Sample = NULL,
                             balance.sample = TRUE,
                             nfeature = 2000,
                             nfeature.per.se = 50,
                             min.cells = 20,
                             downsample = 2500,
                             ncores = 1,
                             seed = 2024){
  if(!(is.matrix(scdata) | is(scdata, "sparseMatrix"))){
    scdata = as.matrix(scdata)
  }
  scmeta$CellType = scmeta[, CellType]
  scmeta$SE = scmeta[, SE]
  if(nrow(scmeta)!=ncol(scdata)) stop("Please ensure that the rows in scmeta matches the columns in scdata.")
  if(min(scdata)>=0 & max(scdata)>80) scdata = log2(scdata+1)

  cts <- scmeta %>% count(CellType, SE) %>% filter(n>min.cells) %>%
    count(CellType) %>% filter(n>1) %>% pull(CellType)
  Ws <- mclapply(cts, function(x){
    tmpmeta = scmeta[scmeta$CellType==x, ]
    tmpdat = scdata[, scmeta$CellType==x]
    if(is.null(Sample)){
      if(nrow(tmpdat)>nfeature){ # Select variable features
        varfeatures = rownames(tmpdat)[order(-apply(tmpdat, 1, var))]
        varfeatures = varfeatures[1:min(nfeature, length(varfeatures))]
        tmpdat = tmpdat[varfeatures, ]
      }
      if(scale) tmpdat <- as.matrix(Seurat::ScaleData(tmpdat, verbose = FALSE))
    }else{
      tmpmeta$Sample = tmpmeta[, Sample]
      samples = table(tmpmeta$Sample)
      samples = names(samples)[samples>min.cells]

      if(nrow(tmpdat)>nfeature){# select variable genes
        ranks = lapply(samples, function(s){
          rank(apply(tmpdat[, tmpmeta$Sample==s, drop=FALSE], 1, var))
        })
        ranks = Reduce("+", ranks)
        varfeatures = rownames(tmpdat)[order(-ranks)]
        varfeatures = varfeatures[1:nfeature]
        tmpdat = tmpdat[varfeatures, ]
      }
      # univariance normalization
      if(scale){
        tmpdat = lapply(samples, function(s){
          ScaleData(tmpdat[, tmpmeta$Sample==s], verbose = FALSE)
        })
        tmpdat = do.call(cbind, tmpdat)
        tmpmeta = tmpmeta[match(colnames(tmpdat), rownames(tmpmeta)), ]
      }

      # balance the cells from multiple samples
      if(balance.sample){
        set.seed(seed)
        balancesize = max(floor(median(table(tmpmeta$Sample))), min.cells)
        tmpdat = lapply(samples, function(s){
          idx = which(tmpmeta$Sample==s)
          if(length(idx)>balancesize) idx = sample(idx, balancesize)
          tmpdat[, idx]
        })
        tmpdat = do.call(cbind, tmpdat)
        tmpmeta = tmpmeta[match(colnames(tmpdat), rownames(tmpmeta)), ]
      }
    }

    if(nrow(tmpmeta)>downsample){
      set.seed(seed)
      idx = sample(nrow(tmpmeta), downsample)
      tmpmeta = tmpmeta[idx, ]
      tmpdat = tmpdat[, idx]
    }

    #### Generate binary H matrix ####
    H = matrix(0, nrow = length(unique(tmpmeta$SE)), ncol = nrow(tmpmeta),
                dimnames = list(row = unique(tmpmeta$SE), col = colnames(tmpdat)))
    idx = cbind(row = match(tmpmeta$SE, rownames(H)), col = 1:nrow(tmpmeta))
    idx = as.array(idx)
    H[idx] = 1
    H = H[rowSums(H)>2, ]
    if(nrow(H)<2) return(NULL)
    W = NMFGenerateW(t(H), tmpdat, scale = FALSE, nfeature = nfeature,
                     nfeature.per.se = nfeature.per.se)
    W
  }, mc.cores = ncores)
  names(Ws) <- cts
  return(Ws)
}


#' Train SE Deconvolution Model
#'
#' This function trains a Non-negative Matrix Factorization (NMF) model for SE deconvolution
#' based on given spatial ecotype fractions and gene expression matrix. Prior to NMF, each gene
#' is scaled to mean 0 and unit variance. To satisfy the non-negativity requirement of NMF,
#' the expression matrix is processed using posneg transformation, which converts the expression
#' matrix into two matrices, one containing only positive values and the other containing only
#' negative values with the sign inverted. The two matrices are subsequently concatenated to
#' produce the training data.
#'
#' @param Fracs A fraction matrix, with rows as samples.
#' @param ExpMat A gene expression matrix with genes in rows and samples in columns.
#' @param scale Logical indicating whether to scale the gene expression matrix. Default is TRUE.
#' @param nfeature Integer specifying the top variable features for training the models (default: 2000).
#' @param nfeature.per.se Integer specifying the maximal number of features to select for each SE (default: 50).
#' @param method A character string specifying the NMF method to use. Default is "brunet".
#'
#' @return A matrix containing the basis NMF W matrix with rows as features and columns as SEs.
#'
#' @import NMF
#' @export
#'
NMFGenerateW <- function(Fracs, ExpMat, scale = TRUE,
                         nfeature = 2000,
                         nfeature.per.se = 50,
                         method = "brunet"){ # , "nsNMF", "snmf/r", "offset"
  require("NMF")
  if(min(ExpMat)>=0 & max(ExpMat)>80) ExpMat = log2(ExpMat+1)
  if(nrow(ExpMat)>nfeature){
    varfeatures = rownames(ExpMat)[order(-apply(ExpMat, 1, var))]
    varfeatures = varfeatures[1:min(nfeature, length(varfeatures))]
    ExpMat = ExpMat[varfeatures, ]
  }

  if(scale){
    ExpMat <- Seurat::ScaleData(ExpMat, verbose = FALSE)
    ExpMat[is.na(ExpMat)] = 0
    to_predict = ExpMat
    if(sum(to_predict<0)>0) to_predict = posneg(to_predict)
    to_predict[is.na(to_predict)] = 0
  }else{
    to_predict <- as.matrix(ExpMat)
    if(sum(to_predict<0)>0) to_predict = posneg(to_predict)
    to_predict[is.na(to_predict)] = 0
  }
  # NEW CODE ADDED TO FILTER BASED ON VARIANCE
  to_predict = to_predict[apply(to_predict, 1, function(x) var(x) > 0), ]

  Fracs = as.matrix(t(Fracs))
  FracsF = Fracs[,colnames(to_predict)]

  my_method <- function (i, v, x, copy = FALSE, eps = .Machine$double.eps, ...) {
    w <- .basis(x)
    h <- .coef(x)
    nb <- nbterms(x)
    nc <- ncterms(x)
    w <- NMF:::std.divergence.update.w(v, w, h, nbterms = nb, ncterms = nc, copy = copy)
    if (i%%10 == 0) {
      w <- pmax.inplace(w, eps, ibterms(x))
    }
    if (copy) {
      .basis(x) <- w
    }
    return(x)
  }

  dummy = NMF::rnmf(nrow(to_predict), H = FracsF)

  dummyW = dummy@W
  dummyH = dummy@H

  dummyW[is.na(dummyW)] = 0
  dummyH[is.na(dummyH)] = 0

  dummyWF = dummyW
  dummyHF = dummyH[, colnames(to_predict)]


  # NEW CODE ADDED - TOGGLE THIS LINE AS NEEDED
  dummyHF = dummyHF + .Machine$double.eps

  my.seeding.method <- function(model, target){
    basis(model) <- dummyWF #estim.r@fit@W
    # initialize H randomly
    coef(model) <- dummyHF
    # return updated object
    return(model)
  }

  nmf_method <- NMF::NMFStrategy('my-method', method, Update = my_method,
                                 objective = 'KL', Stop='connectivity')

  new_nmf = NMF::nmf(to_predict, nrow(FracsF), nrun = 1, method = nmf_method,
                     seed = my.seeding.method, .opt='P1')

  W = new_nmf@fit@W

  W[is.na(W)] = 0
  colnames(W) = rownames(FracsF)

  idx <- duplicated(rownames(W))
  if(sum(idx)>0){
    rownames(W)[!idx] <- paste0(rownames(W)[!idx], "__pos")
    rownames(W)[idx] <- paste0(rownames(W)[idx], "__neg")
  }

  ## Feature selection
  delta = W + apply(W, 1, function(x){ sort(-x)[2] })
  genes = lapply(1:ncol(delta), function(j){
    x = delta[, j]
    names(x) = rownames(delta)
    x = sort(x, decreasing = TRUE)
    # x = x[grepl("__pos", names(x))]
    genes = gsub("__.*", "", names(x[x>0]))
    if(length(genes)>nfeature.per.se) return(genes[1:nfeature.per.se])
    return(genes)
  })
  names(genes) = colnames(delta)
  idx = gsub("__.*", "", rownames(W)) %in% unlist(genes)
  W = W[idx, ]
  W = W[order(-rowSums(W)), ]
  W
}

