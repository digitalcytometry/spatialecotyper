#' @import NMF
.nmf.predict <- function(W, testdat, scale = FALSE, normalize = FALSE){
  set.seed(39)
  #### to compare across cell type, do not normalize
  testdat <- as.matrix(testdat)
  W <- as.matrix(W)
  if(grepl("__", rownames(W))[1]){
    training_gene_set = intersect(unique(gsub("__.*", "", rownames(W))), rownames(testdat))
  }else{
    training_gene_set = intersect(rownames(W), rownames(testdat))
  }
  if(length(training_gene_set)<5){
    H0 = matrix(NA, nrow = ncol(testdat), ncol = ncol(W))
    rownames(H0) = colnames(testdat)
    colnames(H0) = colnames(W)
    warning("Too few genes (<5) overlap with the model. Return NA predictions. ")
    return(t(H0))
  }
  testdat = testdat[match(training_gene_set, rownames(testdat)), , drop = FALSE]
  if(scale){ testdat <- Seurat::ScaleData(testdat, verbose = FALSE) }
  testdat[is.na(testdat)] = 0
  rownames(testdat) = training_gene_set
  testdat <- as.matrix(testdat)
  if(grepl("__", rownames(W))[1]){
    to_predict = as.matrix(NMF::posneg(testdat))
    idx <- duplicated(rownames(to_predict))
    rownames(to_predict)[!idx] <- paste0(rownames(to_predict)[!idx], "__pos")
    rownames(to_predict)[idx] <- paste0(rownames(to_predict)[idx], "__neg")
  }else{
    to_predict = as.matrix(testdat)
  }
  to_predict = to_predict[apply(to_predict, 1, function(x) var(x) > 0), , drop = FALSE]
  features <- intersect(rownames(to_predict), rownames(W))
  W <- W[match(features, rownames(W)), , drop = FALSE]
  to_predict <- to_predict[match(features, rownames(to_predict)), , drop = FALSE]
  W[is.na(W)] <- 0
  to_predict[is.na(to_predict)] <- 0

  my_method <- function (i, v, x, copy = FALSE, eps = .Machine$double.eps, ...){
    w <- .basis(x)
    h <- .coef(x)
    nb <- nbterms(x)
    nc <- ncterms(x)
    h <- NMF:::std.divergence.update.h(v, w, h, nbterms = nb, ncterms = nc, copy = copy)
    #w <- NMF:::std.divergence.update.w(v, w, h, nbterms = nb, ncterms = nc, copy = copy)
    if (i%%10 == 0) {
      h <- pmax.inplace(h, eps, icterms(x))
      #w <- pmax.inplace(w, eps, ibterms(x))
    }
    if (copy) {
      #.basis(x) <- w
      .coef(x) <- h
    }
    return(x)
  }

  ws = as.matrix(W)
  dummy = NMF::rnmf(ncol(ws), to_predict)

  my.seeding.method <- function(model, target){
    basis(model) <- ws #estim.r@fit@W
    # initialize H randomly
    coef(model) <- dummy@H
    # return updated object
    return(model)
  }

  nmf_method <- NMF::NMFStrategy('my-method', 'brunet', Update = my_method, objective = 'KL', Stop='connectivity')
  new_nmf = NMF::nmf(to_predict, ncol(ws), nrun = 1, method = nmf_method, seed = my.seeding.method, .opt='P1')

  H = new_nmf@fit@H
  rownames(H) = colnames(ws)
  if(normalize) H <- t(t(H) / colSums(H))
  return(H)
}

#' Prediction Using Pretrained NMF Model
#'
#' This function uses pretrained NMF models to recover cell states / spatial ecotypes.
#' It takes a factorization matrix W representing a pretrained NMF model and a numeric gene expression matrix.
#'
#' @param W Matrix representing the factorization matrix W of a pretrained NMF model.
#' @param testdat Numeric matrix containing the new data for which NMF scores are to be predicted.
#' @param scale Logical indicating whether to scale the input data.
#' @param ncell.per.run Integer specifying the maximum number of cells per NMF prediction run to avoid memory issues.
#' @param sum2one Logical. If `TRUE`, normalizes the predicted scores so that
#' they sum to 1 for each sample.
#' @param ncores Integer specifying the number of CPU cores to use for parallel processing.
#'
#' @return A matrix representing the NMF prediction scores, with rows as samples/cells and columns as SE groups.
#'
#' @examples
#' library(data.table)
#' library(SpatialEcoTyper)
#' bulkdata <- fread("https://spatialecotyper.stanford.edu/inc/inc.public.vignettes.php?file=SKCM_RNASeqV2.geneExp.tsv",
#'                   sep = "\t", header = TRUE, data.table = FALSE)
#' rownames(bulkdata) = bulkdata[, 1]
#' bulkdata = as.matrix(bulkdata[, -1])
#' W <- readRDS(file.path(system.file("extdata", package = "SpatialEcoTyper"), "Bulk_SE_Recovery_W.rds"))
#'
#' # Predict SE abundances in bulk tumors
#' preds <- NMFpredict(W = W, bulkdata, scale = TRUE)
#' head(preds[, 1:5])
#'
#' @importFrom parallel mclapply
#' @export

NMFpredict <- function(W, testdat,
                       scale = FALSE,
                       ncell.per.run = 500,
                       sum2one = TRUE,
                       ncores = 1){

  ## Prediction
  if(scale){
    ngenes = colSums(testdat>1e-16)
    cantpred = colnames(testdat)[ngenes<3]
    if(length(cantpred)>0){
      warning(length(cantpred), " cells are omitted due to a lack of genes.")
    }
    testdat = testdat[, ngenes>=3, drop = FALSE]
    testdat = Seurat::ScaleData(testdat, verbose = FALSE)
  }
  if(ncol(testdat)>ncell.per.run){
    # Split cells into as-equal-as-possible contiguous chunks (sizes differ by at most 1).
    nfold <- ceiling(ncol(testdat)/ncell.per.run)
    base_size <- ncol(testdat) %/% nfold
    remainder <- ncol(testdat) %% nfold
    fold_sizes <- rep(base_size, nfold)
    if(remainder>0) fold_sizes[seq_len(remainder)] <- fold_sizes[seq_len(remainder)] + 1
    fold_ends <- cumsum(fold_sizes)
    fold_starts <- fold_ends - fold_sizes + 1
    H <- mclapply(1:nfold, function(x){
      idx <- fold_starts[x]:fold_ends[x]
      tmpdat2 <- testdat[, idx, drop = FALSE]
      .nmf.predict(W, tmpdat2, scale = FALSE, normalize = sum2one)
    }, mc.cores = ncores)
    H <- do.call(cbind, H)
  }else{
    H = .nmf.predict(W, testdat, scale = FALSE, normalize = sum2one)
  }
  return(t(H))
}
