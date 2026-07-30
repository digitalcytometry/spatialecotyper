#' Select the Optimal Rank from a Cophenetic Coefficient Profile
#'
#' Selects the best rank (e.g., number of NMF clusters/communities) from a
#' series of candidate ranks based on their cophenetic coefficients. Among
#' ranks whose cophenetic coefficient exceeds a stability threshold, this
#' returns the rank followed by the largest and most persistent subsequent
#' decline in cophenetic coefficient, indicating a transition to less
#' stable solutions at higher ranks.
#'
#' @param cophenetic Numeric vector of cophenetic coefficients, one per rank
#' in \code{ranks}, ordered so that \code{cophenetic[i]} corresponds to
#' \code{ranks[i]}.
#' @param ranks Numeric vector of candidate ranks corresponding to
#' \code{cophenetic}. Defaults to \code{seq_along(cophenetic)} (i.e., ranks
#' are assumed to be given in order starting from 1). Internally sorted in
#' increasing order together with \code{cophenetic} before use, so callers
#' do not need to pre-sort.
#' @param threshold Numeric cophenetic coefficient cutoff (default: 0.95).
#' Only ranks with a cophenetic coefficient strictly greater than this value
#' are considered as candidates.
#' @param min_rank Numeric minimum rank value to consider as a candidate
#' (default: 2), regardless of its cophenetic coefficient. Use this to rule
#' out trivially small ranks, which can otherwise score well simply because
#' cophenetic coefficients tend to be close to 1 at very low ranks.
#' @param persistence_window Integer number of ranks ahead over which to
#' measure the \emph{persistent} decline in cophenetic coefficient
#' (default: 3). For a candidate rank \code{i}, this is the average
#' per-rank drop from rank \code{i} to rank \code{i + persistence_window}
#' -- i.e. \code{(cophenetic[i] - cophenetic[i + persistence_window]) /
#' persistence_window} -- which distinguishes a decline that is sustained
#' over several ranks from a single-step wobble that may recover
#' immediately afterward.
#' @param immediate_weight Numeric weight in \code{[0, 1]} (default: 0.5)
#' balancing the immediate next-rank drop against the persistence-window
#' drop when scoring candidates (see Details). A value of 1 uses only the
#' immediate next-rank drop; a value of 0 uses only the persistence-window
#' drop.
#'
#' @return The selected best rank (a single value from \code{ranks}): the
#' candidate rank (cophenetic coefficient above \code{threshold} and rank
#' at least \code{min_rank}) with the highest combined drop score (see
#' Details). If no rank satisfies both conditions, a warning is issued and
#' the rank with the overall highest cophenetic coefficient is returned
#' instead.
#'
#' @details
#' For each rank \code{i}, two measures of subsequent decline are computed:
#' \describe{
#'   \item{\code{immediate_drop}}{The drop in cophenetic coefficient from
#'   rank \code{i} to the very next rank, \code{cophenetic[i] -
#'   cophenetic[i + 1]}.}
#'   \item{\code{persistence_drop}}{The average per-rank drop over the next
#'   \code{persistence_window} ranks, \code{(cophenetic[i] - cophenetic[i +
#'   persistence_window]) / persistence_window}. Near the end of the
#'   profile, where \code{i + persistence_window} exceeds the number of
#'   ranks tested, this falls back to \code{immediate_drop} instead.}
#' }
#' These are combined into a single score for each candidate,
#' \code{immediate_weight * immediate_drop + (1 - immediate_weight) *
#' persistence_drop}, and the candidate with the highest score is returned.
#' Unlike a purely peak-based approach, a candidate rank does not need to be
#' a strict local maximum of the cophenetic profile -- any rank above
#' \code{threshold} and \code{min_rank} is scored on how much, and how
#' persistently, the cophenetic coefficient falls apart afterward, so this
#' also works when the coefficient declines steadily without oscillating.
#'
#' Because both \code{immediate_drop} and \code{persistence_drop} are
#' undefined (\code{NA}) for the very last rank in \code{ranks} (there is no
#' subsequent rank to measure a decline against), its score is always
#' \code{NA} and it will only be returned if it happens to be the sole
#' candidate.
#'
#' @export
SelectBestRank <- function(cophenetic,
                           ranks = seq_along(cophenetic),
                           threshold = 0.95,
                           min_rank = 2,
                           persistence_window = 3,
                           immediate_weight = 0.5) {

  stopifnot(length(ranks) == length(cophenetic))
  cophenetic <- cophenetic[order(ranks)]
  ranks <- sort(ranks)
  n <- length(cophenetic)
  if (n == 1) return(ranks[1])

  # Drop in cophenetic coefficient from rank i to the next rank i + 1.
  # The last rank has no "next" rank, so its drop is undefined (NA).
  immediate_drop <- c(cophenetic[-n] - cophenetic[-1], NA)
  persistence_drop <- rep(NA_real_, n)

  for (i in seq_len(n - 1)) {
    persistence_drop[i] <- (
      cophenetic[i] - cophenetic[i + persistence_window]
    ) / persistence_window
  }
  persistence_drop[is.na(persistence_drop)] =immediate_drop[is.na(persistence_drop)]
  candidate <- which(cophenetic > threshold & ranks >= min_rank)

  if (length(candidate) == 0) {
    warning("No rank satisfies the cophenetic coefficient threshold.")
    return(ranks[which.max(cophenetic)])
  }

  score <- (immediate_weight * immediate_drop[candidate] +
              (1 - immediate_weight) * persistence_drop[candidate])

  # Among qualifying candidates, prefer the highest rank.
  result <- data.frame(rank = ranks[candidate], cophenetic = cophenetic[candidate],
                       score = score)
  result <- result[order(-result$score), ]
  return(result$rank[1])
}
