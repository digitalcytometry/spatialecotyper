#' Assess concordance between true and predicted spatial ecotype labels
#'
#' Computes an F1-score matrix comparing true spatial ecotype (SE) labels
#' and predicted SE labels (e.g., from cross-validation). The function
#' calculates precision and recall per SE within each sample and returns
#' an averaged F1-score matrix across samples.
#'
#' @param scmeta A data.frame containing single-cell metadata.
#' @param SE Character. Column name in `scmeta` for true SE labels.
#' @param Pred Character. Column name in `scmeta` for predicted SE labels.
#' @param CellType Character or NULL. Optional column name for cell type.
#' If provided, SE labels will be concatenated with cell type labels.
#' @param Sample Character or NULL. Column name for sample identifiers in `scmeta`.
#' If NULL, all cells are treated as coming from a single sample.
#' @param metric One of 'F1', 'F2', 'precision', or 'recall'.
#'
#' @return A matrix of averaged F1 scores with rows and columns corresponding
#' to SE labels.
#'
#' @details
#' The function computes the specified metrics (e.g. F1 scores) for all SEs,
#' and then average them across samples.
#'
#' @importFrom dplyr group_by count mutate arrange
#' @importFrom tidyr pivot_wider
#' @importFrom rlang .data
#' @export
#'
ComputeMetrics <- function(scmeta, SE = "SE",
                           Pred = "cvPred",
                           CellType = NULL,
                           Sample = "Sample",
                           metric = c("F1", "F2", "precision", "recall")[1]){
  if(is.null(Sample)){
    scmeta$Sample = "Sample"
  }else{
    scmeta$Sample = scmeta[, Sample]
  }
  if(!is.null(CellType)){
    scmeta = scmeta %>% arrange(.data[[CellType]], .data[[SE]])
    scmeta[, SE] = paste0(scmeta[, SE], "_", scmeta[, CellType])
    scmeta[, Pred] = paste0(scmeta[, Pred], "_", scmeta[, CellType])
  }
  ses = unique(scmeta[, SE])
  metric_list = lapply(unique(scmeta$Sample), function(ss){
    metas = scmeta[scmeta$Sample==ss, ]
    Recalls <- metas %>% group_by(.data[[SE]]) %>% count(.data[[Pred]]) %>%
      mutate(Fracs = n / sum(n)) %>%
      pivot_wider(id_cols = all_of(SE), names_from = all_of(Pred), values_from = Fracs)
    Recalls <- as.data.frame(Recalls)
    rownames(Recalls) <- Recalls[,1]
    Recalls <- Recalls[, -1]
    Recalls <- as.matrix(Recalls)
    Recalls <- Recalls[match(ses, rownames(Recalls)), match(ses, colnames(Recalls))]
    rownames(Recalls) <- ses
    colnames(Recalls) <- ses
    Recalls[is.na(Recalls)] = 0
    Recalls <- Recalls / rowSums(Recalls)

    Precision <- metas %>% group_by(.data[[Pred]]) %>% count(.data[[SE]]) %>%
      mutate(Fracs = n / sum(n)) %>%
      pivot_wider(id_cols = all_of(Pred), names_from = all_of(SE), values_from = Fracs) %>%
      as.data.frame()
    rownames(Precision) <- Precision[, 1]
    Precision <- as.matrix(Precision[, -1])
    Precision <- Precision[match(ses, rownames(Precision)), match(ses, colnames(Precision))]
    rownames(Precision) <- ses
    colnames(Precision) <- ses
    Precision[is.na(Precision)] <- 0
    Precision = Precision / rowSums(Precision)

    # ---- Metric selection ----
    if (metric == "recall") {
      out <- Recalls
    } else if (metric == "precision") {
      out <- Precision
    } else {
      Recalls[is.na(Recalls)] = 0
      Precision[is.na(Precision)] <- 0
      beta <- ifelse(metric == "F2", 2, 1)
      out <- (1 + beta^2) * Precision * Recalls /
        (beta^2 * Precision + Recalls)
    }
    return(out)
  })
  result = Reduce("+", lapply(metric_list, function(x){ x[is.na(x)] = 0; x })) /
    Reduce("+", lapply(metric_list, function(x){ !is.na(x) }))
  result[is.na(result)] = 0
  return(result)
}
