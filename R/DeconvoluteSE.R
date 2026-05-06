#' Infer SE Abundances Using a Pretrained NMF Model
#'
#' This function estimates the abundances of SEs in bulk tumor gene expression data
#' using a provided basis matrix (`W`) from NMF model.
#'
#' @param dat A numeric matrix of gene expression data, e.g. TPM matrix for bulk RNA-seq data, logCPM matrix for Visium data etc.
#' @param scale Logical. If `TRUE`, the input data is scaled before making predictions.
#' @param W A numeric matrix representing the basis matrix (`W`) from a pretrained NMF model.
#' @param nsample.per.run Integer specifying the maximum number of samples to process per
#' prediction run to manage memory usage efficiently.
#' @param sum2one Logical. If `TRUE`, normalizes the predicted SE abundances so that
#' they sum to 1 for each sample.
#' @param ncores Integer specifying the number of CPU cores to use for parallel processing.
#'
#' @return A matrix of SE abundances in bulk tumors.
#'
#' @examples
#'
#' bulkdata <- fread("https://spatialecotyper.stanford.edu/inc/inc.public.vignettes.php?file=SKCM_RNASeqV2.geneExp.tsv",
#'                   sep = "\t", header = TRUE, data.table = FALSE)
#' rownames(bulkdata) = bulkdata[, 1]
#' bulkdata = as.matrix(bulkdata[, -1])
#'
#' # Predict SE abundances in bulk tumors
#' se_abundances <- DeconvoluteSE(dat = bulkdata)
#' head(se_abundances[, 1:5])
#'
#' @export
#'
#' @importFrom parallel mclapply
#'
DeconvoluteSE <- function(dat, scale = TRUE, W = NULL,
                          nsample.per.run = 500,
                          sum2one = TRUE, ncores = 8){
  if(is.null(W)){
    W <- readRDS(file.path(system.file("extdata", package = "SpatialEcoTyper"),
                           "Bulk_SE_Recovery_W.rds"))
  }
  genes = unique(gsub("_.*", "", rownames(W)))
  dat = dat[rownames(dat) %in% genes, ]
  if(all(dat>=0) & max(dat)>80){
    dat = log2(dat+1)
  }
  dat = dat[rownames(dat) %in% gsub("_.*", "", rownames(W)), ]
  if(scale) dat = ScaleData(dat, verbose = FALSE)
  # if(min(dat)>=0){
  #   ngenes = colSums(dat<1e-16)
  #   cantpred = colnames(dat)[ngenes<3]
  #   if(length(cantpred)>0){
  #     warning(length(cantpred), " samples are omitted due to lack of model gene expression")
  #   }
  # }
  preds <- NMFpredict(W = W, dat, scale = FALSE,
                      ncell.per.run = nsample.per.run,
                      sum2one = sum2one, ncores = ncores)
  return(preds)
}
