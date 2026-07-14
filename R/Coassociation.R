#' Compute co-association of cell states across samples
#'
#' This function computes pairwise co-association (Pearson correlation)
#' between cell states (defined as combinations of spatial ecotype (SE) and
#' cell type) across samples. Co-association is evaluated under multiple
#' abundance calculation schemes that differ in the inclusion of NonSE cells
#' and the treatment of absent cell states, and the resulting correlation
#' matrices are averaged to obtain a robust co-association matrix.
#'
#' @param scmeta A data.frame containing single-cell metadata.
#' @param Sample Character. Column name in \code{scmeta} specifying sample IDs.
#' @param SE Character. Column name in \code{scmeta} specifying spatial ecotype labels.
#' @param CellType Character. Column name in \code{scmeta} specifying cell type annotations.
#' @param NonSE Character. Label used to denote non-SE cells (default: \code{"NonSE"}).
#' @param nperm Integer. Number of permutations used for significance testing
#'   when \code{test = TRUE}. Default is 1000.
#' @param test Logical. If \code{TRUE} (default), performs permutation testing
#'   using \code{CoassociationTest()} and returns both the co-association matrix
#'   and the corresponding p-value matrix. If \code{FALSE}, only the
#'   co-association matrix is returned.
#'
#' @return
#' If \code{test = TRUE}, a list with the following components:
#' \describe{
#'   \item{CoassociationIndex}{A symmetric matrix of averaged Pearson
#'   correlation coefficients between cell states.}
#'   \item{Pval}{A matrix of permutation-based p-values for the
#'   co-association scores.}
#' }
#'
#' If \code{test = FALSE}, a symmetric matrix of averaged Pearson correlation
#' coefficients between cell states.
#'
#' @details
#' Cell states are defined as concatenations of SE and cell type labels.
#' For each sample and cell type, the relative abundance of each cell state
#' is calculated under four schemes:
#' \enumerate{
#'   \item Including both SE and NonSE cell states, with absent states represented as \code{NA}.
#'   \item Excluding NonSE cell states, with absent states represented as \code{NA}.
#'   \item Including both SE and NonSE cell states, with absent states filled with zero.
#'   \item Excluding NonSE cell states, with absent states filled with zero.
#' }
#'
#' Pairwise Pearson correlations between cell states are then computed across
#' samples for each abundance matrix using
#' \code{cor(method = "pearson", use = "pairwise.complete.obs")}. The final
#' co-association matrix is obtained by averaging the correlation coefficients
#' across the four schemes while ignoring missing values.
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
