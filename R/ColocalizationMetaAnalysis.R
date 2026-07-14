#' Meta-analysis of colocalization results across samples
#'
#' Combines per-sample pairwise colocalization index matrices and per-SE
#' p-values into a single meta-analysis summary, using Stouffer's method for
#' both the colocalization index matrix and the SE-level statistics.
#'
#' @param colocalization_results list, one element per sample. Each element
#'   must itself be a list containing:
#'   \describe{
#'     \item{ColocIndex}{a matrix of pairwise colocalization indices,
#'       with `rownames()`/`colnames()` equal to cell state labels.}
#'     \item{Pval}{a named numeric vector of per-SE two-sided p-values, as
#'       returned by \code{\link{Colocalization}}/\code{\link{CoassociationTest}}.
#'       When these p-values carry a `"Zscore"` attribute (the default when
#'       produced by those functions), the signed Z-scores are used to
#'       combine samples so that effect direction is preserved; otherwise the
#'       Z-score magnitude is recovered from the p-value and direction is
#'       treated as unknown.}
#'   }
#' @param cap numeric > 0, symmetric cap applied to `ColocIndex` values
#'   before combining (values are clipped to `[-cap, cap]`) to limit the
#'   influence of extreme outlier indices from any single sample.
#' @param min.samples integer >= 1, minimum number of samples a cell state
#'   must be backed by to be retained in the combined result.
#'   Cell states appearing in the `ColocIndex` matrix of fewer than
#'   `min.samples` samples are dropped.
#'
#' @return A list with:
#'   \describe{
#'     \item{MetaColocIndex}{a states x states matrix of integrated
#'       colocalization indices.}
#'     \item{MetaPval}{a named numeric vector of combined per-SE
#'       two-sided p-values, Stouffer-combined across samples.}
#'   }
#'
#' @importFrom dplyr group_by summarize mutate
#' @importFrom stats qnorm pnorm
#' @export
ColocalizationMetaAnalysis <- function(colocalization_results, cap = 5,
                                       min.samples = 1) {

  if (!is.list(colocalization_results) || length(colocalization_results) == 0) {
    stop("`colocalization_results` must be a non-empty list.")
  }
  if (!(cap > 0)) stop("`cap` must be a positive number.")
  if (!(min.samples >= 1)) stop("`min.samples` must be >= 1.")

  has_fields <- vapply(colocalization_results, function(xx) {
    !is.null(xx$ColocIndex) && !is.null(xx$Pval)
  }, logical(1))
  if (!all(has_fields)) {
    stop("Every element of `colocalization_results` must contain both ",
         "`ColocIndex` and `Pval`.")
  }

  ## ---- 1. Combine pairwise ColocIndex matrices ----
  # cap extreme values per sample
  mat_list <- lapply(colocalization_results, function(xx) {
    tmp <- as.matrix(xx$ColocIndex)
    tmp[tmp >  cap] <-  cap
    tmp[tmp < -cap] <- -cap
    tmp
  })

  # keep states observed in at least `min.samples` samples
  state_count <- table(unlist(lapply(mat_list, rownames)))
  states <- names(state_count[state_count >= min.samples])
  if (length(states) < 3) {
    stop("Fewer than three cell states are shared across samples ",
         "(after applying `min.samples`).")
  }

  # reorder/pad every sample's matrix to the common `states` set
  # (states absent from a given sample become NA rows/columns for that
  # sample only, and are excluded cell-by-cell below)
  mat_list <- lapply(mat_list, function(x) {
    x <- x[match(states, rownames(x)), match(states, colnames(x))]
    dimnames(x) <- list(states, states)
    x
  })

  mat_array <- simplify2array(mat_list)
  metagg <- apply(mat_array, c(1, 2), sum, na.rm = TRUE)

  na_array     <- simplify2array(lapply(mat_list, is.na))
  n_nonmissing <- length(mat_list) - apply(na_array, c(1, 2), sum)

  metagg <- metagg / (sqrt(n_nonmissing) + 1e-7)
  metagg[n_nonmissing == 0] <- NA  # no contributing samples -> unknown, not 0

  ## ---- 2. Combine per-state p-values (Stouffer's method) ----

  eps <- .Machine$double.eps
  pval_df <- lapply(colocalization_results, function(xx) {
    zscore <- attr(xx$Pval, "Zscore")
    if (is.null(zscore)) {
      p <- pmin(pmax(as.numeric(xx$Pval), eps), 1 - eps)
      zscore <- qnorm(p / 2, lower.tail = FALSE)
    }
    data.frame(SE = names(xx$Pval), Zscore = as.numeric(zscore))
  })
  pval_df <- do.call(rbind, pval_df)
  pval_df <- pval_df[!is.na(pval_df$Zscore), ]

  meta_pval <- pval_df %>%
    dplyr::group_by(SE) %>%
    dplyr::summarize(n = length(Zscore),
                     Zscore = sum(Zscore) / (sqrt(n) + 1e-7)) %>%
    dplyr::mutate(Pval = pnorm(abs(Zscore), lower.tail = FALSE) * 2) %>%
    as.data.frame()
  Pval <- meta_pval$Pval
  names(Pval) <- meta_pval$SE
  list(MetaColocIndex = metagg, MetaPval = Pval)
}
