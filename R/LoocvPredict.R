#' Perform leave-one-out cross-validation (LOOCV) for SE prediction
#'
#' Trains NMF-based spatial ecotype (SE) recovery models using subsets of
#' single-cell data and evaluates predictions on held-out data. Supports
#' repeated cross-validation and parallel execution.
#'
#' @param scdata A gene expression matrix (genes x cells).
#' @param scmeta A data.frame containing metadata for each cell.
#' Row names of `scmeta` should match the column names in `scdata`.
#' @param Sample Character. Column name in `scmeta` specifying sample IDs.
#' If fewer than two unique samples are present, cells are randomly split
#' into training and test sets within each cell type.
#' @param CellType Character. Column name specifying cell type annotations in `scmeta`.
#' @param SE Character. Column name specifying spatial ecotype labels in `scmeta`.
#' @param repeats Integer. Number of cross-validation repeats.
#' @param ncores Integer. Number of cores for parallel computation.
#' @param scale Boolean specifying whether to perform univariance normalization
#' for training and validation data (default: TRUE).
#' @param verbose Boolean specifying whether to print the log messages.
#' @param ... Additional arguments passed to \code{NMFGenerateWList}.
#'
#' @return The input `scmeta` data.frame with an added column `cvPred`
#' containing predicted SE labels for each cell.
#'
#' @details
#' For each repeat:
#' \itemize{
#'   \item If multiple samples are available, performs leave-one-sample-out CV
#'   \item Otherwise, performs random stratified splitting within cell types
#'   \item Trains NMF models using \code{NMFGenerateWList}
#'   \item Predicts SE labels using \code{RecoverSE}
#' }
#' Predictions across repeats are aggregated by majority vote per cell.
#'
#' @importFrom parallel mclapply
#' @importFrom dplyr group_by mutate ungroup count arrange distinct
#' @importFrom rlang .data
#' @export
#'
LoocvPredict = function(scdata, scmeta,
                        Sample = "Sample",
                        CellType = "CellType",
                        SE = "SE",
                        repeats = 30,
                        ncores = 4,
                        scale = TRUE,
                        verbose = TRUE,
                        ...){
  seeds = sample(1:10000, repeats)

  preds = parallel::mclapply(1:repeats, function(ii){
    set.seed(seeds[ii])

    if(length(unique(scmeta[[Sample]]))<2){
      scmeta = scmeta %>% dplyr::group_by(.data[[CellType]]) %>%
        mutate(Split = sample(c(rep("train", floor(n()/2)),
                                rep("test", n() - floor(n()/2))))) %>%
      ungroup()
    }else{
      scmeta$Split = scmeta[[Sample]]
    }

    preds = lapply(unique(scmeta$Split), function(ss){
      message(Sys.time(), " Training on ", ss)
      trainMeta = scmeta[scmeta$Split==ss, ]
      trainDat = scdata[, scmeta$Split==ss]
      testMeta = scmeta[scmeta$Split!=ss, ]
      testDat = scdata[, scmeta$Split!=ss]
      if(verbose){
        Ws = NMFGenerateWList(trainDat, trainMeta,
                              CellType = CellType,
                              SE = SE, scale = scale,
                              Sample = "Split",
                              seed = seeds[ii], ...)
        preds = RecoverSE(dat = testDat,
                          celltypes = testMeta$CellType,
                          ncell.per.run = 5000,
                          Ws = Ws, min.score = 0,
                          scale = scale)
      }else{
        suppressMessages(Ws = NMFGenerateWList(trainDat, trainMeta,
                                               CellType = CellType,
                                               SE = SE, scale = scale,
                                               Sample = "Split",
                                               seed = seeds[ii], ...))
        suppressMessages(preds = RecoverSE(dat = testDat,
                                           celltypes = testMeta$CellType,
                                           ncell.per.run = 5000,
                                           Ws = Ws, min.score = 0,
                                           scale = scale))
      }
      preds
    })
    preds = do.call(rbind, preds)
    preds
  }, mc.cores = ncores)
  preds = do.call(rbind, preds)
  preds = preds %>% count(CID, SE) %>% arrange(desc(n)) %>%
    distinct(CID, .keep_all = TRUE)
  scmeta$cvPred = preds$SE[match(colnames(scdata), preds$CID)]
  scmeta
  return(scmeta)
}


