#' Aggregate recovery models across runs
#'
#' Aggregate multiple recovery model weight matrices by identifying
#' consistently selected marker genes across models and averaging
#' their weights.
#'
#' This function is designed for aggregating recovery models generated
#' from repeated training runs or cross-validation folds. For each cell
#' type, genes are retained if they are consistently selected across
#' models based on:
#'
#' \itemize{
#'   \item a minimum delta threshold required for a gene to be selected
#'   \item a minimum fraction of models in which the gene is selected
#' }
#'
#' Positive and negative feature pairs (`__pos` and `__neg`) are retained
#' for all selected genes. Final weights are computed as the mean across
#' non-missing values.
#'
#' @param model_list A list of recovery model objects. Each element should contain a
#' named list of weight matrices indexed by cell type, derived from `NMFGenerateWList`.
#'
#' @param delta.threshold Minimum delta weights required for a gene to be
#' considered associated with a spatial ecotype. Default is `0.01`.
#'
#' @param min.model.fraction Minimum fraction of models in which a gene
#' must be selected to be retained. Default is `0.5`.
#'
#' @return A named list of aggregated weight matrices, one per cell type.
#'
#' @examples
#' \dontrun{
#' Ws_list = lapply(1:30, function(ii){
#'   Ws = NMFGenerateWList(scdata, scmeta,
#'                         CellType = "CellType",
#'                         SE = "SE", Sample = "Sample",
#'                         nfeature = 300,
#'                         nfeature.per.se = 50,
#'                         ncores = 8)
#' })
#' aggregated_models <- AggregateRecoverModels(
#'   model_list = Ws_list,
#'   delta.threshold = 0.01,
#'   min.model.fraction = 0.5
#' )
#' }
#'
#' @export
#'
AggregateRecoverModels = function(model_list,
                                  delta.threshold = 0.01,
                                  min.model.fraction = 0.5){
  Ws_final = lapply(names(model_list[[1]]), function(ct){
    Ws = lapply(model_list, function(xx){ xx[[ct]] })
    gene_df = lapply(Ws, function(W){
      delta = W + apply(W, 1, function(x){ sort(-x)[2] })
      genes <- lapply(1:ncol(delta), function(j){
        x <- delta[, j]
        names(x) <- rownames(delta)
        x <- sort(x[x>delta.threshold], decreasing = TRUE)
        x <- x[grepl("__pos", names(x))]
        gs <- gsub("_.*", "", names(x))
        return(gs)
      })
      data.frame(SE = rep(colnames(delta), lengths(genes)),
                 Gene = unlist(genes))
    })
    gene_df = do.call(rbind, gene_df)
    gene_df = gene_df %>% count(SE, Gene) %>%
      mutate(frac = n / length(model_list)) %>% filter(frac>min.model.fraction)
    features = unique(c(paste0(gene_df$Gene, "__pos"), paste0(gene_df$Gene, "__neg")))
    ses = unique(unlist(lapply(Ws, colnames)))
    Ws = lapply(Ws, function(xx){
      xx = xx[match(features, rownames(xx)), match(ses, colnames(xx))]
      xx
    })
    W <- Reduce("+", lapply(Ws, function(x){ x[is.na(x)]=0; x})) /
      Reduce("+", lapply(Ws, function(x){ !is.na(x) }))
    rownames(W) <- features
    colnames(W) <- ses
    W
  })
  names(Ws_final) = names(model_list[[1]])
  return(Ws_final)
}
