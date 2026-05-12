#' Infer Cell Numbers per Spot from Expression Profiles
#'
#' Estimates the number of cells for each Visium spot based on
#' total gene expression. The inferred values are linearly scaled such
#' that the spot with the lowest total expression is assigned 1 cell,
#' and the average inferred number of cells across all spots is approximately
#' equal to `avg.number`.
#'
#' @param normdat A numeric matrix or data frame with genes in rows and
#'   spots in columns. Values should represent normalized or raw
#'   expression counts.
#' @param avg.number Numeric value specifying the desired average number
#'   of inferred cells per spot. Default is 5 for Visium.
#'
#' @return An integer vector containing the inferred number of cells for
#'   each spot
#'
#' @export
#'
InferNCells <- function(normdat, avg.number = 5) {

  if (!is.numeric(avg.number) || length(avg.number) != 1 || avg.number <= 0) {
    stop("'avg.number' must be a single positive numeric value.")
  }
  total <- colSums(normdat, na.rm = TRUE)

  if (length(total) == 0) {
    stop("Input contains no columns.")
  }

  if (all(total == total[1])) {
    return(rep(round(avg.number), length(total)))
  }

  slope <- (avg.number - 1) / (mean(total) - min(total))
  intercept <- 1 - slope * min(total)

  ncells <- round(slope * total + intercept)

  # Ensure at least 1 cell per spot
  ncells[ncells < 1] <- 1

  return(as.integer(ncells))
}
