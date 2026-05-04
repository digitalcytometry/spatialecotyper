#' Compute co-association of cell states across samples
#'
#' This function computes pairwise co-association (Pearson correlation)
#' between cell states (defined as combinations of SE and cell type)
#' across samples. It evaluates co-association under multiple data
#' transformations (with/without non-SE states and missing value handling)
#' and returns an averaged correlation matrix.
#'
#' @param scmeta A data.frame containing single-cell metadata.
#' @param Sample Character. Column name in \code{scmeta} specifying sample IDs.
#' @param SE Character. Column name in \code{scmeta} specifying spatial ecotype labels.
#' @param CellType Character. Column name in \code{scmeta} specifying cell type annotations.
#' @param NonSE Character. Label used to denote non-SE cells (default: "NonSE").
#'
#' @return A symmetric matrix of Pearson correlation coefficients between
#' cell states, averaged across multiple normalization strategies.
#'
#' @details
#' Cell states are defined as concatenations of SE and cell type labels.
#' The function computes abundances of SE-associated cell state under four schemes:
#' (1) including all SE and NonSE states, (2) excluding NonSE states,
#' (3) same as in (1), but treating zero abundance as missing values (NA), and
#' (4) same as in (2), but treating zero abundance as NA. Pairwise Pearson correlations
#' between cell states were then calculated across all scRNA-seq samples for each abundance
#' matrix, using the cor function in R with pairwise complete observations. The final co-association
#' values were obtained by averaging the correlations across the four schemes.
#'
#' @importFrom dplyr count group_by mutate filter
#' @importFrom tidyr pivot_wider
#' @importFrom rlang .data
#' @export
Coassociation = function(scmeta, Sample = "Sample", SE = "SE", CellType = "CellType",
                         NonSE = "NonSE", nperm = 1000, test = TRUE){
  # ---- Define cell states ----
  scmeta$CellState <- paste0(scmeta[[SE]], "_", scmeta[[CellType]])
  states <- sort(unique(scmeta$CellState))

  # ---- Helper function to compute fraction matrix ----
  get_state_matrix <- function(meta, drop_nonSE = FALSE, fill_zero = FALSE) {
    if (drop_nonSE) { meta <- meta[meta[[SE]] != NonSE, ] }
    df <- meta %>% dplyr::count(.data[[Sample]], .data[[CellType]], .data$CellState) %>%
      dplyr::group_by(.data[[Sample]], .data[[CellType]]) %>%
      dplyr::mutate(Frac = n / sum(n)) %>%
      tidyr::pivot_wider(id_cols = all_of(Sample), names_from = "CellState",
                         values_from = "Frac")
    df <- as.data.frame(df)
    rownames(df) <- df[[Sample]]
    mat <- as.matrix(df[, -1, drop = FALSE])
    mat <- mat[, match(states, colnames(mat)), drop = FALSE]
    if (fill_zero) { mat[is.na(mat)] <- 0 }
    return(mat)
  }

  # ---- Generate matrices under different settings ----
  F1 <- get_state_matrix(scmeta, drop_nonSE = FALSE, fill_zero = FALSE)
  F2 <- get_state_matrix(scmeta, drop_nonSE = TRUE,  fill_zero = FALSE)
  F3 <- get_state_matrix(scmeta, drop_nonSE = FALSE, fill_zero = TRUE)
  F4 <- get_state_matrix(scmeta, drop_nonSE = TRUE,  fill_zero = TRUE)

  F_list <- list(F1, F2, F3, F4)

  # ---- Compute correlations ----
  cor_list <- suppressWarnings(
    lapply(F_list, function(mat) {
      cor(mat, method = "pearson", use = "pairwise.complete.obs")
    })
  )

  # ---- Average correlations ----
  num <- Reduce("+", lapply(cor_list, function(x) { x[is.na(x)] <- 0; x }))
  denom <- Reduce("+", lapply(cor_list, function(x) { !is.na(x) }))
  denom[denom == 0] <- NA

  scgg <- num / denom
  scgg[is.na(scgg)] <- 0

  if(test){
    Pval = CoassociationTest(scgg, nperm = nperm)
    return(list(CoassociationIndex = scgg, Pval = Pval))
  }else{
    return(scgg)
  }
}
